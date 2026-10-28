#!/usr/bin/env python3
"""Download GFS (NOMADS) or ECMWF IFS open data and write it in grib2json format for earth.

Usage:
  scripts/fetch-weather.py gfs 925 surface                # latest analysis -> data/weather/current/
  scripts/fetch-weather.py ecmwf 925 surface --timeline   # also past 48 h + 72 h forecast, 3-hourly
  scripts/fetch-weather.py gfs 925 --cycle 2026102800     # pin the "latest" cycle

Writes wind (U/V) and relative humidity for each level, plus mean sea level pressure and 3-hour
precipitation (mm), all on a 1.0 deg grid. Levels are isobaric levels in hPa, or "surface"
(10 m wind, 2 m relative humidity).
Requires: pip install eccodes
"""
import argparse
import collections
import concurrent.futures
import datetime as dt
import json
import os
import sys
import time
import urllib.error
import urllib.request

import eccodes
import numpy as np

WEATHER_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "public", "data", "weather"))
GFS_URL = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl"
ECMWF_URLS = ["https://data.ecmwf.int/forecasts", "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com"]
CENTER_NAMES = {7: "US National Weather Service - NCEP(WMC)", 98: "European Centre for Medium-Range Weather Forecasts"}
STEP = dt.timedelta(hours=3)
CYCLE = dt.timedelta(hours=6)

Field = collections.namedtuple("Field", "short level start header values")


def http_get(url, headers=None):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers or {}), timeout=120) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return None
        except Exception:
            pass
        time.sleep(2 * (attempt + 1))
    return None


def gfs_download(cycle, fhour, levels):
    lev = ["lev_%s_mb=on" % l for l in levels if l != "surface"]
    if "surface" in levels:
        lev += ["lev_10_m_above_ground=on", "lev_2_m_above_ground=on"]
    lev += ["lev_mean_sea_level=on", "lev_surface=on"]
    url = "%s?file=gfs.t%02dz.pgrb2.1p00.f%03d&%s&var_UGRD=on&var_VGRD=on&var_RH=on&var_PRMSL=on&var_APCP=on" \
          "&dir=%%2Fgfs.%s%%2F%02d%%2Fatmos" % (GFS_URL, cycle.hour, fhour, "&".join(lev),
                                               cycle.strftime("%Y%m%d"), cycle.hour)
    data = http_get(url)
    return data if data and data.startswith(b"GRIB") else None


def ecmwf_download(cycle, fhour, levels):
    stream = "oper" if cycle.hour in (0, 12) else "scda"
    path = "%s/%02dz/ifs/0p25/%s/%s-%dh-%s-fc" % (
        cycle.strftime("%Y%m%d"), cycle.hour, stream, cycle.strftime("%Y%m%d%H0000"), fhour, stream)
    isobaric = {str(l) for l in levels if l != "surface"}
    surface = {"msl", "tp"} | ({"10u", "10v", "2t", "2d"} if "surface" in levels else set())
    for base in ECMWF_URLS:
        index = http_get("%s/%s.index" % (base, path))
        if not index:
            continue
        wanted = []
        for line in index.decode().splitlines():
            r = json.loads(line)
            if (r["levtype"] == "pl" and r.get("levelist") in isobaric and r["param"] in ("u", "v", "r")) or \
                    (r["levtype"] == "sfc" and r["param"] in surface):
                wanted.append(r)
        ranges = [{"Range": "bytes=%d-%d" % (r["_offset"], r["_offset"] + r["_length"] - 1)} for r in wanted]
        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            chunks = list(pool.map(lambda h: http_get("%s/%s.grib2" % (base, path), h), ranges))
        if chunks and all(chunks):
            return b"".join(chunks)
    return None


def decode(grib):
    """Returns a Field for each GRIB2 message, subsampled to a 1.0 deg grid."""
    fields, offset = [], 0
    while offset < len(grib):
        length = int.from_bytes(grib[offset + 8:offset + 16], "big")
        gid = eccodes.codes_new_from_message(grib[offset:offset + length])
        offset += length
        get = lambda key, ktype=None: eccodes.codes_get(gid, key, ktype=ktype) if ktype else eccodes.codes_get(gid, key)
        increment = get("iDirectionIncrementInDegrees")
        k = max(1, int(round(1.0 / increment)))
        values = eccodes.codes_get_values(gid).reshape(get("Nj"), get("Ni"))[::k, ::k]
        center = get("centre", int)
        header = {
            "discipline": get("discipline"), "parameterCategory": get("parameterCategory"),
            "parameterNumber": get("parameterNumber"), "parameterUnit": get("units"),
            "surface1Type": get("typeOfFirstFixedSurface", int), "surface1Value": get("level"),
            "forecastTime": get("endStep", int), "scanMode": 0, "center": center,
            "centerName": CENTER_NAMES.get(center, get("centre")),
            "refTime": "%sT%02d:00:00.000Z" % (
                dt.datetime.strptime(str(get("dataDate")), "%Y%m%d").strftime("%Y-%m-%d"), get("dataTime") // 100),
            "nx": values.shape[1], "ny": values.shape[0], "numberPoints": values.size,
            "lo1": get("longitudeOfFirstGridPointInDegrees"), "la1": get("latitudeOfFirstGridPointInDegrees"),
            "dx": increment * k, "dy": get("jDirectionIncrementInDegrees") * k,
        }
        header["lo2"] = (header["lo1"] + (header["nx"] - 1) * header["dx"]) % 360
        header["la2"] = header["la1"] - (header["ny"] - 1) * header["dy"]
        if header["surface1Type"] == 100:
            header["surface1Value"] *= 100  # grib2json reports isobaric levels in Pa
        fields.append(Field(get("shortName"), get("level"), get("startStep", int), header, values))
        eccodes.codes_release(gid)
    return fields


def record(field, decimals, name, values=None, **header):
    values = field.values if values is None else values
    data = np.round(values.ravel().astype(float), decimals)
    return {"header": dict(field.header, parameterNumberName=name, **header),
            "data": (data.astype(int) if decimals == 0 else data).tolist()}


def to_products(fields, levels):
    """Maps decoded fields to {(type, layer): [records]} for wind, relative humidity and MSLP."""
    by = {}
    for f in fields:
        by.setdefault((f.short, f.level), f)
    out = {}
    for level in levels:
        if level == "surface":
            layer, u, v, rh = "surface-level", by.get(("10u", 10)), by.get(("10v", 10)), by.get(("2r", 2))
            t, td = by.get(("2t", 2)), by.get(("2d", 2))
            if rh is None and t is not None and td is not None:  # ECMWF has no 2 m RH: derive it (Magnus formula)
                tc, tdc = t.values - 273.15, td.values - 273.15
                values = np.clip(100 * np.exp(17.625 * tdc / (243.04 + tdc) - 17.625 * tc / (243.04 + tc)), 0, 100)
                out[("relative_humidity", layer)] = [record(t, 0, "Relative_humidity", values,
                                                            parameterCategory=1, parameterNumber=1, parameterUnit="%")]
        else:
            layer = "isobaric-%shPa" % level
            u, v, rh = by.get(("u", int(level))), by.get(("v", int(level))), by.get(("r", int(level)))
        if u is not None and v is not None:
            out[("wind", layer)] = [record(u, 2, "U-component_of_wind"), record(v, 2, "V-component_of_wind")]
        if rh is not None:
            out[("relative_humidity", layer)] = [record(rh, 0, "Relative_humidity")]
    msl = by.get(("prmsl", 0)) or by.get(("msl", 0))
    if msl is not None:
        out[("mean_sea_level_pressure", None)] = [record(msl, 0, "Pressure_reduced_to_MSL")]
    return out


class Source:
    def __init__(self, model, levels):
        self.model, self.levels = model, levels
        self.download = gfs_download if model == "gfs" else ecmwf_download
        self.cache = collections.OrderedDict()

    def get(self, cycle, fhour):
        """Returns the decoded fields for a cycle and forecast hour, or None if unavailable."""
        key = (cycle, fhour)
        if key not in self.cache:
            grib = self.download(cycle, fhour, self.levels)
            self.cache[key] = decode(grib) if grib else None
            while len(self.cache) > 8:
                self.cache.popitem(last=False)
            time.sleep(0.3)  # be gentle with rate limits
        return self.cache[key]

    def precipitation(self, cycle, fhour):
        """3-hour precipitation (mm) ending at cycle+fhour, from GFS buckets or ECMWF since-start totals."""
        tp = lambda fields: next((f for f in fields or [] if f.short == "tp"), None)
        end = tp(self.get(cycle, fhour))
        if end is None:
            return None
        scale = 1000.0 if end.header["parameterUnit"] == "m" else 1.0
        values = end.values * scale
        if end.start != fhour - 3:
            prev = tp(self.get(cycle, fhour - 3))
            if prev is None or prev.start != end.start:
                return None
            values = values - prev.values * scale
        return record(end, 2, "Total_precipitation_3h", np.clip(values, 0, None),
                      parameterUnit="kg m-2", forecastTime=fhour)


def write(products, model, valid=None):
    """Writes products to data/weather/current/ (valid=None) or data/weather/yyyy/mm/dd/hhmm-...json."""
    for (kind, layer), records in products.items():
        stamp, folder = ("current", "current") if valid is None else (valid.strftime("%H%M"), valid.strftime("%Y/%m/%d"))
        path = os.path.join(WEATHER_DIR, folder, "-".join(filter(None, [stamp, kind, layer, model, "1.0"])) + ".json")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(records, f, separators=(",", ":"))


def fetch_slot(src, valid, latest, dated):
    """Writes all products valid at `valid` from the newest cycle that has them. Returns False if none found."""
    cycle = min(latest, valid.replace(hour=valid.hour // 6 * 6))
    for _ in range(2):  # fall back to the previous cycle if a file is not yet published
        fhour = int((valid - cycle).total_seconds() // 3600)
        fields = src.get(cycle, fhour)
        if fields:
            break
        cycle -= CYCLE
    else:
        return False
    products = to_products(fields, src.levels)
    precip = src.precipitation(cycle, fhour) if fhour >= 3 else None
    for back in (1, 2):  # analysis steps have no precipitation: use an earlier cycle's forecast instead
        if precip is None and fhour < 3:
            precip = src.precipitation(cycle - back * CYCLE, fhour + back * 6)
    if precip:
        products[("precipitation", None)] = [precip]
    write(products, src.model, valid if dated else None)
    return True


def latest_cycle(src):
    cycle = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0, tzinfo=None)
    cycle = cycle.replace(hour=cycle.hour // 6 * 6)
    for _ in range(10):
        if src.get(cycle, 0):
            return cycle
        cycle -= CYCLE
    raise RuntimeError("no %s cycle found in the last 60 hours" % src.model)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("model", choices=["gfs", "ecmwf"])
    p.add_argument("levels", nargs="+", help='isobaric levels in hPa (e.g. 925) and/or "surface"')
    p.add_argument("--cycle", help="YYYYMMDDHH to use as the latest cycle (default: newest available)")
    p.add_argument("--timeline", action="store_true", help="also fetch past and forecast hours (3-hourly)")
    p.add_argument("--past-hours", type=int, default=48)
    p.add_argument("--forecast-hours", type=int, default=72)
    a = p.parse_args()
    src = Source(a.model, a.levels)
    latest = dt.datetime.strptime(a.cycle, "%Y%m%d%H") if a.cycle else latest_cycle(src)
    if not fetch_slot(src, latest, latest, dated=False):
        raise RuntimeError("no %s data for %s" % (a.model, latest))
    print("wrote current %s %s (%s)" % (a.model, ", ".join(a.levels), latest.strftime("%Y-%m-%d %HZ")))
    if a.timeline:
        valid, end, written = latest - dt.timedelta(hours=a.past_hours), latest + dt.timedelta(hours=a.forecast_hours), 0
        while valid <= end:
            if fetch_slot(src, valid, latest, dated=True):
                written += 1
            else:
                print("  missing %s" % valid.strftime("%Y-%m-%d %H:%MZ"))
            valid += STEP
        print("wrote %d timeline steps (%s to %s)" % (
            written, (latest - dt.timedelta(hours=a.past_hours)).strftime("%Y-%m-%d %HZ"), end.strftime("%Y-%m-%d %HZ")))


if __name__ == "__main__":
    sys.exit(main())
