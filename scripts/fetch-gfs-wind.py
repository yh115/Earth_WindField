#!/usr/bin/env python3
"""Download GFS 1.0 deg U/V wind from NOMADS and write it in grib2json format for earth.

Usage: scripts/fetch-gfs-wind.py 925 [--date YYYYMMDD --cycle HH]
Level is an isobaric level in hPa, or "surface" for 10 m above ground.
Requires: pip install eccodes
"""
import argparse
import datetime as dt
import json
import os
import sys
import urllib.request

import eccodes

FILTER_URL = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl"
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "public", "data", "weather", "current")
HEADER_KEYS = {
    "discipline": "discipline", "center": "centre", "subcenter": "subCentre",
    "parameterCategory": "parameterCategory", "parameterNumber": "parameterNumber",
    "parameterUnit": "units", "forecastTime": "forecastTime",
    "surface1Type": "typeOfFirstFixedSurface", "surface1Value": "level",
    "numberPoints": "numberOfPoints", "scanMode": "scanningMode",
    "nx": "Ni", "ny": "Nj", "lo1": "longitudeOfFirstGridPointInDegrees",
    "la1": "latitudeOfFirstGridPointInDegrees", "lo2": "longitudeOfLastGridPointInDegrees",
    "la2": "latitudeOfLastGridPointInDegrees", "dx": "iDirectionIncrementInDegrees",
    "dy": "jDirectionIncrementInDegrees",
}


def download(date, cycle, level):
    lev = "lev_10_m_above_ground=on" if level == "surface" else "lev_%s_mb=on" % level
    url = "%s?file=gfs.t%sz.pgrb2.1p00.f000&%s&var_UGRD=on&var_VGRD=on&dir=%%2Fgfs.%s%%2F%s%%2Fatmos" % (
        FILTER_URL, cycle, lev, date, cycle)
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = resp.read()
    if not data.startswith(b"GRIB"):
        raise RuntimeError("no GRIB data for %s %sZ" % (date, cycle))
    return data


def latest(level):
    now = dt.datetime.now(dt.timezone.utc)
    for back in range(0, 48, 6):
        t = now - dt.timedelta(hours=back)
        date, cycle = t.strftime("%Y%m%d"), "%02d" % (t.hour // 6 * 6)
        try:
            return date, cycle, download(date, cycle, level)
        except Exception:
            continue
    raise RuntimeError("no GFS cycle found in the last 48 hours")


def to_records(grib):
    records = []
    offset = 0
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
        records.append({"header": header, "data": data})
    records.sort(key=lambda r: r["header"]["parameterNumber"])  # U (2) before V (3)
    return records


def main():
    p = argparse.ArgumentParser()
    p.add_argument("level", help='isobaric level in hPa (e.g. 925) or "surface"')
    p.add_argument("--date", help="YYYYMMDD (default: latest available)")
    p.add_argument("--cycle", help="00, 06, 12 or 18")
    a = p.parse_args()
    if a.date and a.cycle:
        date, cycle, grib = a.date, a.cycle, download(a.date, a.cycle, a.level)
    else:
        date, cycle, grib = latest(a.level)
    name = "surface-level" if a.level == "surface" else "isobaric-%shPa" % a.level
    out = os.path.join(OUT_DIR, "current-wind-%s-gfs-1.0.json" % name)
    with open(out, "w") as f:
        json.dump(to_records(grib), f, separators=(",", ":"))
    print("wrote %s (GFS %s %sZ)" % (os.path.normpath(out), date, cycle))


if __name__ == "__main__":
    sys.exit(main())
