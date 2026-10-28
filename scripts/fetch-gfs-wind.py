#!/usr/bin/env python3
"""Download GFS 1.0 deg U/V wind from NOMADS and write it in grib2json format for earth.

Usage:
  scripts/fetch-gfs-wind.py 925 surface                 # latest analysis -> data/weather/current/
  scripts/fetch-gfs-wind.py 925 surface --timeline      # also past 48 h + 72 h forecast, 3-hourly
  scripts/fetch-gfs-wind.py 925 --date 20261028 --cycle 00

Levels are isobaric levels in hPa, or "surface" for 10 m above ground.
Requires: pip install eccodes
"""
import argparse
import datetime as dt
import json
import os
import sys
import time
import urllib.request

import eccodes

FILTER_URL = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl"
WEATHER_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "public", "data", "weather"))
HEADER_KEYS = {
    "discipline": "discipline", "subcenter": "subCentre",
    "parameterCategory": "parameterCategory", "parameterNumber": "parameterNumber",
    "parameterUnit": "units", "forecastTime": "forecastTime", "surface1Value": "level",
    "numberPoints": "numberOfPoints", "scanMode": "scanningMode",
    "nx": "Ni", "ny": "Nj", "lo1": "longitudeOfFirstGridPointInDegrees",
    "la1": "latitudeOfFirstGridPointInDegrees", "lo2": "longitudeOfLastGridPointInDegrees",
    "la2": "latitudeOfLastGridPointInDegrees", "dx": "iDirectionIncrementInDegrees",
    "dy": "jDirectionIncrementInDegrees",
}
STEP = dt.timedelta(hours=3)


def download(cycle, fhour, levels):
    """Returns GRIB2 bytes for the given cycle (datetime) and forecast hour, or None if unavailable."""
    lev = "&".join("lev_10_m_above_ground=on" if l == "surface" else "lev_%s_mb=on" % l for l in levels)
    url = "%s?file=gfs.t%02dz.pgrb2.1p00.f%03d&%s&var_UGRD=on&var_VGRD=on&dir=%%2Fgfs.%s%%2F%02d%%2Fatmos" % (
        FILTER_URL, cycle.hour, fhour, lev, cycle.strftime("%Y%m%d"), cycle.hour)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                data = resp.read()
            return data if data.startswith(b"GRIB") else None
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return None
        except Exception:
            pass
        time.sleep(2 * (attempt + 1))
    return None


def latest_cycle(levels):
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0, tzinfo=None)
    cycle = now.replace(hour=now.hour // 6 * 6)
    for _ in range(8):
        grib = download(cycle, 0, levels)
        if grib:
            return cycle, grib
        cycle -= dt.timedelta(hours=6)
    raise RuntimeError("no GFS cycle found in the last 48 hours")


def layer_name(header):
    return "surface-level" if header["surface1Type"] == 103 else "isobaric-%dhPa" % (header["surface1Value"] // 100)


def to_layers(grib):
    """Splits GRIB2 bytes into {layer name: [U record, V record]} in grib2json format."""
    layers, offset = {}, 0
    while offset < len(grib):
        length = int.from_bytes(grib[offset + 8:offset + 16], "big")
        gid = eccodes.codes_new_from_message(grib[offset:offset + length])
        offset += length
        header = {k: eccodes.codes_get(gid, v) for k, v in HEADER_KEYS.items()}
        header["center"] = eccodes.codes_get(gid, "centre", ktype=int)
        header["surface1Type"] = eccodes.codes_get(gid, "typeOfFirstFixedSurface", ktype=int)
        if header["surface1Type"] == 100:
            header["surface1Value"] *= 100  # grib2json reports isobaric levels in Pa
        header["refTime"] = "%sT%02d:00:00.000Z" % (
            dt.datetime.strptime(str(eccodes.codes_get(gid, "dataDate")), "%Y%m%d").strftime("%Y-%m-%d"),
            eccodes.codes_get(gid, "dataTime") // 100)
        header["centerName"] = "US National Weather Service - NCEP(WMC)"
        header["parameterNumberName"] = "U-component_of_wind" if header["parameterNumber"] == 2 else "V-component_of_wind"
        data = [round(float(x), 2) for x in eccodes.codes_get_values(gid)]
        eccodes.codes_release(gid)
        layers.setdefault(layer_name(header), []).append({"header": header, "data": data})
    for records in layers.values():
        records.sort(key=lambda r: r["header"]["parameterNumber"])  # U (2) before V (3)
    return layers


def write_layers(grib, valid=None):
    """Writes each layer to data/weather/current/ (valid=None) or data/weather/yyyy/mm/dd/hhmm-...json."""
    for name, records in to_layers(grib).items():
        if valid is None:
            path = os.path.join(WEATHER_DIR, "current", "current-wind-%s-gfs-1.0.json" % name)
        else:
            path = os.path.join(WEATHER_DIR, valid.strftime("%Y/%m/%d"),
                                "%s-wind-%s-gfs-1.0.json" % (valid.strftime("%H%M"), name))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(records, f, separators=(",", ":"))


def fetch_timeline(latest, levels, past_hours, forecast_hours):
    """Fills 3-hourly slots from latest-past_hours to latest+forecast_hours using the newest cycle for each."""
    valid, end, written = latest - dt.timedelta(hours=past_hours), latest + dt.timedelta(hours=forecast_hours), 0
    while valid <= end:
        cycle = min(latest, valid.replace(hour=valid.hour // 6 * 6))
        for _ in range(2):  # fall back to the previous cycle if a file is not yet published
            fhour = int((valid - cycle).total_seconds() // 3600)
            grib = download(cycle, fhour, levels)
            if grib:
                write_layers(grib, valid)
                written += 1
                break
            cycle -= dt.timedelta(hours=6)
        else:
            print("  missing %s" % valid.strftime("%Y-%m-%d %H:%MZ"))
        time.sleep(0.5)  # be gentle with NOMADS rate limits
        valid += STEP
    return written


def main():
    p = argparse.ArgumentParser()
    p.add_argument("levels", nargs="+", help='isobaric levels in hPa (e.g. 925) and/or "surface"')
    p.add_argument("--date", help="YYYYMMDD (default: latest available)")
    p.add_argument("--cycle", help="00, 06, 12 or 18")
    p.add_argument("--timeline", action="store_true", help="also fetch past and forecast hours (3-hourly)")
    p.add_argument("--past-hours", type=int, default=48)
    p.add_argument("--forecast-hours", type=int, default=72)
    a = p.parse_args()
    if a.date and a.cycle:
        cycle = dt.datetime.strptime(a.date + a.cycle, "%Y%m%d%H")
        grib = download(cycle, 0, a.levels)
        if not grib:
            raise RuntimeError("no GRIB data for %s %sZ" % (a.date, a.cycle))
    else:
        cycle, grib = latest_cycle(a.levels)
    write_layers(grib)
    print("wrote current %s (GFS %s)" % (", ".join(a.levels), cycle.strftime("%Y-%m-%d %HZ")))
    if a.timeline:
        n = fetch_timeline(cycle, a.levels, a.past_hours, a.forecast_hours)
        print("wrote %d timeline steps (%s to %s)" % (
            n, (cycle - dt.timedelta(hours=a.past_hours)).strftime("%Y-%m-%d %HZ"),
            (cycle + dt.timedelta(hours=a.forecast_hours)).strftime("%Y-%m-%d %HZ")))


if __name__ == "__main__":
    sys.exit(main())
