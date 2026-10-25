#!/usr/bin/env python3
"""Build one CARMA alumina event netCDF file per droplet size from the 2025 reentries.

Interpreter: any python3 with numpy + scipy (the netCDF file is written with
scipy.io.netcdf_file, so netCDF4 is not needed).

Inputs:
  - The reentry inventory (2025_emissions.csv), read as in
    reentry_csv_to_alumina_emissions.py: Al2O3 in tonnes, reentry_lat_deg,
    reentry_lon_deg; rows with no Al2O3 or no location are skipped.
  - Droplet Emission Profiles/droplet_<D>um_satellite_al_aerosol_profile.csv
    (from run_droplet_emission_profiles.py):
       Altitude,Surviving Al Mass Fraction
    Fraction of the released mass that remains as aerosol (after
    vaporisation) in each 1 km bin at the droplet's final altitude
    (Altitude is the bin centre), for a single droplet diameter D. It sums
    to the surviving share (about 0.67).

Mass basis:
  The inventory Al2O3 is the total released mass, so each reentry emits
      Al2O3 * droplet(z)
  at altitude z. --conserve-mass instead normalises each droplet profile to 1
  (shape only, vaporisation loss ignored, total mass unchanged).

  Only satellite droplet profiles exist. By default rocket bodies (S0-S4) use
  the satellite profile too; --satellites-only drops them.

Output (same layout as csv_to_alumina_events.py, 64-bit offset netCDF):
  <outdir>/alumina_2025_reentry_events_d<D>um.nc, one row per reentry and
  non-zero profile bin; dates shifted by --year-offset (default -2024 -> year 1).

Usage:
  droplet_profiles_to_alumina_events.py
  droplet_profiles_to_alumina_events.py --conserve-mass -o other_dir
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

from reentry_csv_to_alumina_emissions import (TONNE_TO_KG, InputError, format_date,
                                              profile_for, read_profile)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROFILE_DIR = ROOT / "Droplet Emission Profiles"
PROFILE_RE = re.compile(r"^droplet_([0-9.]+)um_satellite_al_aerosol_profile\.csv$")

# Noleap calendar: February always has 28 days.
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def size_tag(d_um):
    """0.01 -> 'd0p01um', 10 -> 'd10um' (same names as psd_sensitivitytest/)."""
    return "d" + f"{d_um:g}".replace(".", "p") + "um"


def model_date(start, year_offset, where):
    """'YYYY-MM-DD HH:MM:SS' -> (YYYYMMDD, seconds of day) in the noleap calendar."""
    ymd, hms = start.split(" ")
    year, month, day = (int(x) for x in ymd.split("-"))
    h, m, s = (int(x) for x in hms.split(":"))
    year += year_offset
    if year < 0:
        raise InputError(f"{where}: year {year} < 0 after --year-offset {year_offset}")
    if month == 2 and day == 29:
        raise InputError(f"{where}: Feb 29 does not exist in the noleap calendar")
    if not (1 <= month <= 12 and 1 <= day <= DAYS_IN_MONTH[month - 1]):
        raise InputError(f"{where}: invalid date {start!r}")
    return year * 10000 + month * 100 + day, h * 3600 + m * 60 + s


def read_events(path, satellites_only, year_offset):
    """Located reentries with Al2O3 > 0, in chronological order.

    Returns [(start, date, datesec, lat, lon, mass_kg, category)] and a dict of
    skipped counts/masses.
    """
    events, skipped = [], {}

    def skip(key, mass_t):
        n, m = skipped.get(key, (0, 0.0))
        skipped[key] = (n + 1, m + mass_t)

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        needed = ["date", "category", "Al2O3", "reentry_lat_deg", "reentry_lon_deg"]
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            raise InputError(f"{path}: missing column(s) {missing}")
        for lineno, r in enumerate(reader, start=2):
            where = f"{path}:{lineno}"
            cat = r["category"].strip()
            mass_t = float(r["Al2O3"] or 0.0)
            prof = profile_for(cat)
            if prof is None:
                skip(f"unknown category {cat}", mass_t)
                continue
            if mass_t <= 0.0:
                skip("zero Al2O3", mass_t)
                continue
            lat_s, lon_s = r["reentry_lat_deg"].strip(), r["reentry_lon_deg"].strip()
            if not lat_s or not lon_s:
                skip("no location", mass_t)
                continue
            if satellites_only and prof != "satellite":
                skip("rocket body (--satellites-only)", mass_t)
                continue
            lat, lon = float(lat_s), float(lon_s)
            if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 360.0):
                raise InputError(f"{where}: lat/lon {lat}, {lon} out of range")
            start = format_date(r["date"], where)
            date, datesec = model_date(start, year_offset, where)
            events.append((start, date, datesec, lat, lon % 360.0,
                           mass_t * TONNE_TO_KG, prof))
    # Stable sort on the date string, as in reentry_csv_to_alumina_emissions.py,
    # so event ids match alumina_2025_reentry_events.nc.
    events.sort(key=lambda e: e[0])
    return events, skipped


def write_nc(out_path, events, alt, weights, attrs):
    """One row per (event, non-zero bin); variables as in csv_to_alumina_events.py."""
    keep = weights > 0.0
    alt_k, w_k = alt[keep], weights[keep]
    nlev = len(alt_k)
    n = len(events) * nlev

    eid = np.repeat(np.arange(1, len(events) + 1, dtype="i4"), nlev)
    date = np.repeat(np.array([e[1] for e in events], dtype="i4"), nlev)
    sec = np.repeat(np.array([e[2] for e in events], dtype="i4"), nlev)
    lat = np.repeat(np.array([e[3] for e in events], dtype="f8"), nlev)
    lon = np.repeat(np.array([e[4] for e in events], dtype="f8"), nlev)
    altitude = np.tile(alt_k.astype("f8"), len(events))
    mass = np.outer([e[5] for e in events], w_k).ravel()

    with netcdf_file(out_path, "w", version=2) as nc:
        nc.createDimension("row", n)

        def var(name, dtype, units, long_name, values):
            v = nc.createVariable(name, dtype, ("row",))
            if units:
                v.units = units
            v.long_name = long_name
            v[:] = values

        var("event_id", "i4", "", "emission event id", eid)
        var("date", "i4", "YYYYMMDD", "event start date (noleap model calendar)", date)
        var("datesec", "i4", "s", "event start time, seconds of day", sec)
        var("lat", "f8", "degrees_north", "event latitude", lat)
        var("lon", "f8", "degrees_east", "event longitude in [0, 360)", lon)
        var("altitude", "f8", "km", "release altitude above ground", altitude)
        var("mass", "f8", "kg", "emitted Al2O3 mass", mass)

        total = math.fsum(mass)
        # scipy stores bare Python floats as float32 attributes; force double.
        for k, v in attrs.items():
            setattr(nc, k, np.float64(v) if isinstance(v, float) else v)
        nc.total_mass_kg = np.float64(total)
        nc.calendar = "noleap"
        nc.n_events = len(events)
        nc.history = "created by droplet_profiles_to_alumina_events.py"
    return n, total


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reentries", type=Path,
                   default=ROOT / "Historical Reentries/Barker Emission Inventory"
                   / "2025_emissions.csv")
    p.add_argument("--profile-dir", type=Path, default=PROFILE_DIR)
    p.add_argument("-o", "--outdir", type=Path, default=HERE / "droplet_size_events")
    p.add_argument("--year-offset", type=int, default=-2024,
                   help="integer added to every year (default -2024: 2025 -> year 1)")
    p.add_argument("--conserve-mass", action="store_true",
                   help="normalise each droplet profile to 1 instead of keeping "
                        "the surviving share of the released mass")
    p.add_argument("--satellites-only", action="store_true",
                   help="drop rocket bodies (S0-S4) instead of giving them the "
                        "satellite droplet profile")
    args = p.parse_args(argv)

    profiles = []
    for path in args.profile_dir.iterdir():
        m = PROFILE_RE.match(path.name)
        if m:
            profiles.append((float(m.group(1)), path))
    if not profiles:
        raise InputError(f"{args.profile_dir}: no droplet_<D>um_*.csv profiles")
    profiles.sort()

    events, skipped = read_events(args.reentries, args.satellites_only, args.year_offset)
    inventory_kg = math.fsum(e[5] for e in events)
    print(f"{args.reentries.name}: {len(events)} located reentries, "
          f"{inventory_kg:,.2f} kg Al2O3 released")
    for k, (n, m) in sorted(skipped.items()):
        print(f"  skipped {n:5d} ({m * TONNE_TO_KG:12,.2f} kg): {k}")

    mass_basis = ("droplet profile normalised to 1 (inventory Al2O3 conserved)"
                  if args.conserve_mass else
                  "inventory Al2O3 x droplet profile (fraction of released mass "
                  "remaining as aerosol after vaporisation)")
    rocket_bodies = ("excluded" if args.satellites_only else
                     "included, using the satellite droplet profile")

    args.outdir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'D_um':>6} {'rows':>7} {'surviving':>9} {'total_kg':>14}  file")
    for d_um, path in profiles:
        alt, frac = read_profile(path)
        if len(np.unique(alt)) != len(alt):
            raise InputError(f"{path}: repeated altitudes")
        surviving = math.fsum(frac)
        if not 0.0 < surviving <= 1.0 + 1e-9:
            raise InputError(f"{path}: profile sums to {surviving}; expected the "
                             "surviving share of the released mass, in (0, 1]")
        weights = frac / surviving if args.conserve_mass else frac

        attrs = {
            "title": (f"CARMA alumina point emission events, 2025 reentries, "
                      f"assumed particle diameter {d_um:g} um"),
            "particle_diameter_um": d_um,
            "source_inventory": str(args.reentries),
            "ablation_profile": str(path),
            "mass_basis": mass_basis,
            "rocket_bodies": rocket_bodies,
            "surviving_share": surviving,
            "altitude_note": "altitude is the centre of a 1 km profile bin; "
                             "bins with zero mass are omitted",
            "year_offset": args.year_offset,
        }
        out = args.outdir / f"alumina_2025_reentry_events_{size_tag(d_um)}.nc"
        n, total = write_nc(out, events, alt, weights, attrs)
        print(f"{d_um:6g} {n:7d} {surviving:9.4f} {total:14,.2f}  {out.name}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except InputError as e:
        sys.exit(f"error: {e}")
