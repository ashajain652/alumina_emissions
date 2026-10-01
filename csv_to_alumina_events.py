#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Convert an alumina emissions CSV into the CARMA alumina event netCDF file.

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
(numpy + netCDF4; the miniforge base env has neither).

Input CSV (see ALUMINA_CARMA_PLAN.md section 4):
  id,start_date,lat,long,altitude_km,emitted_mass_kg
  1,0001-01-01 02:24:00,16.1918,-170.9961,85,5.1745

Output netCDF (classic 64-bit offset format, for Fortran/PIO readers):
  dim row
  event_id(row) int, date(row) int YYYYMMDD, datesec(row) int,
  lat(row), lon(row) [0, 360), altitude(row) km above ground, mass(row) kg
  global attrs: source_csv, total_mass_kg, calendar = "noleap"

Dates are parsed by hand (no datetime), because year 1 and the noleap
calendar are not handled consistently by datetime/strftime.

Usage:
  csv_to_alumina_events.py alumina_test_emissions.csv alumina_test_events.nc
  csv_to_alumina_events.py in.csv out.nc --year-offset -2023 --alt-min 20 --alt-max 120
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# Noleap calendar: February always has 28 days.
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

DATE_RE = re.compile(
    r"^\s*(-?\d{1,6})-(\d{1,2})-(\d{1,2})"          # YYYY-MM-DD
    r"(?:[ T](\d{1,2}):(\d{2})(?::(\d{2}(?:\.\d*)?))?)?\s*$"  # [ HH:MM[:SS[.fff]]]
)

REQUIRED_COLUMNS = ["id", "start_date", "lat", "long", "altitude_km", "emitted_mass_kg"]


class InputError(Exception):
    pass


def parse_date(text, year_offset, where):
    """Parse 'YYYY-MM-DD HH:MM:SS' in the noleap calendar.

    Returns (date_int YYYYMMDD, datesec). Rejects Feb 29 and invalid fields.
    """
    m = DATE_RE.match(text)
    if not m:
        raise InputError(f"{where}: cannot parse start_date {text!r} "
                         "(expected YYYY-MM-DD HH:MM:SS)")
    year = int(m.group(1)) + year_offset
    month = int(m.group(2))
    day = int(m.group(3))
    hour = int(m.group(4) or 0)
    minute = int(m.group(5) or 0)
    sec_f = float(m.group(6) or 0.0)

    if year < 0:
        raise InputError(f"{where}: year {year} < 0 after --year-offset {year_offset}")
    if not 1 <= month <= 12:
        raise InputError(f"{where}: month {month} out of range in {text!r}")
    if month == 2 and day == 29:
        raise InputError(f"{where}: Feb 29 in {text!r} does not exist in the noleap "
                         "calendar; shift the date or use --year-offset")
    if not 1 <= day <= DAYS_IN_MONTH[month - 1]:
        raise InputError(f"{where}: day {day} out of range for month {month} in {text!r}")
    if not (0 <= hour <= 23 and 0 <= minute <= 59 and 0.0 <= sec_f < 60.0):
        raise InputError(f"{where}: invalid time of day in {text!r}")
    if sec_f != int(sec_f):
        raise InputError(f"{where}: fractional seconds are not supported ({text!r})")

    date_int = year * 10000 + month * 100 + day
    datesec = hour * 3600 + minute * 60 + int(sec_f)
    return date_int, datesec


def read_csv(path, year_offset, alt_min, alt_max):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        header = [h.strip() for h in (reader.fieldnames or [])]
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise InputError(f"{path}: missing column(s) {missing}; found {header}")
        for lineno, raw in enumerate(reader, start=2):
            r = {k.strip(): (v.strip() if v is not None else "") for k, v in raw.items()}
            if not any(r.values()):
                continue  # blank line
            where = f"{path}:{lineno}"
            try:
                eid = int(r["id"])
                lat = float(r["lat"])
                lon = float(r["long"])
                alt = float(r["altitude_km"])
                mass = float(r["emitted_mass_kg"])
            except ValueError as e:
                raise InputError(f"{where}: bad number ({e})")
            date_int, datesec = parse_date(r["start_date"], year_offset, where)

            for name, v in (("lat", lat), ("long", lon), ("altitude_km", alt),
                            ("emitted_mass_kg", mass)):
                if not math.isfinite(v):
                    raise InputError(f"{where}: {name} is not finite")
            if not -90.0 <= lat <= 90.0:
                raise InputError(f"{where}: lat {lat} outside [-90, 90]")
            if not -180.0 <= lon <= 360.0:
                raise InputError(f"{where}: long {lon} outside [-180, 360]")
            if not alt_min <= alt <= alt_max:
                raise InputError(f"{where}: altitude {alt} km outside "
                                 f"[{alt_min}, {alt_max}] (--alt-min/--alt-max)")
            if mass < 0.0:
                raise InputError(f"{where}: negative mass {mass}")

            lon360 = lon % 360.0
            rows.append(dict(line=lineno, id=eid, date=date_int, datesec=datesec,
                             lat=lat, lon=lon360, alt=alt, mass=mass,
                             mass_text=r["emitted_mass_kg"]))
    if not rows:
        raise InputError(f"{path}: no data rows")
    return rows


def validate(rows, path):
    # Chronological order (non-decreasing in time).
    for prev, cur in zip(rows, rows[1:]):
        if (cur["date"], cur["datesec"]) < (prev["date"], prev["datesec"]):
            raise InputError(f"{path}:{cur['line']}: rows not in chronological order "
                             f"({cur['date']:08d} {cur['datesec']}s after "
                             f"{prev['date']:08d} {prev['datesec']}s)")

    # Consistent date/location within each id; ids must be contiguous blocks.
    first = {}
    last_id = None
    for r in rows:
        eid = r["id"]
        if eid in first and eid != last_id:
            raise InputError(f"{path}:{r['line']}: rows of event id {eid} are not "
                             "contiguous")
        if eid not in first:
            first[eid] = r
        else:
            f0 = first[eid]
            for key in ("date", "datesec", "lat", "lon"):
                if r[key] != f0[key]:
                    raise InputError(f"{path}:{r['line']}: event id {eid} has "
                                     f"inconsistent {key} ({r[key]} vs {f0[key]} "
                                     f"on line {f0['line']})")
        last_id = eid

    # Duplicate altitudes within an event are allowed (they simply add up), but warn.
    seen = set()
    for r in rows:
        key = (r["id"], r["alt"])
        if key in seen:
            print(f"WARNING: event {r['id']} has more than one row at {r['alt']} km "
                  "(masses will add)", file=sys.stderr)
        seen.add(key)


def write_nc(rows, out_path, source_csv, total):
    n = len(rows)
    with Dataset(out_path, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.createDimension("row", n)

        def var(name, dtype, units, long_name, values):
            v = nc.createVariable(name, dtype, ("row",))
            if units:
                v.units = units
            v.long_name = long_name
            v[:] = np.asarray(values, dtype=dtype)

        var("event_id", "i4", "", "emission event id", [r["id"] for r in rows])
        var("date", "i4", "YYYYMMDD", "event start date (noleap model calendar)",
            [r["date"] for r in rows])
        var("datesec", "i4", "s", "event start time, seconds of day",
            [r["datesec"] for r in rows])
        var("lat", "f8", "degrees_north", "event latitude", [r["lat"] for r in rows])
        var("lon", "f8", "degrees_east", "event longitude in [0, 360)",
            [r["lon"] for r in rows])
        var("altitude", "f8", "km", "release altitude above ground",
            [r["alt"] for r in rows])
        var("mass", "f8", "kg", "emitted Al2O3 mass", [r["mass"] for r in rows])

        nc.title = "CARMA alumina point emission events"
        nc.source_csv = str(source_csv)
        nc.total_mass_kg = float(total)
        nc.calendar = "noleap"
        nc.n_events = len({r["id"] for r in rows})
        nc.history = "created by csv_to_alumina_events.py"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="input emissions CSV")
    p.add_argument("out", help="output event netCDF")
    p.add_argument("--year-offset", type=int, default=0,
                   help="integer added to every year (map real dates to model years)")
    p.add_argument("--alt-min", type=float, default=0.0,
                   help="minimum allowed altitude, km above ground (default 0)")
    p.add_argument("--alt-max", type=float, default=140.0,
                   help="maximum allowed altitude, km above ground (default 140)")
    args = p.parse_args(argv)

    src = Path(args.csv).resolve()
    try:
        rows = read_csv(src, args.year_offset, args.alt_min, args.alt_max)
        validate(rows, src)
    except InputError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    total = math.fsum(r["mass"] for r in rows)
    write_nc(rows, args.out, src, total)

    n_events = len({r["id"] for r in rows})
    print(f"Wrote {args.out}: {len(rows)} rows, {n_events} events")
    print(f"  first {rows[0]['date']:08d} {rows[0]['datesec']}s, "
          f"last {rows[-1]['date']:08d} {rows[-1]['datesec']}s (noleap)")
    print(f"  altitude range {min(r['alt'] for r in rows)}-"
          f"{max(r['alt'] for r in rows)} km")
    print(f"  total mass = {total:.6f} kg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
