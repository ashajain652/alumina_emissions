#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Turn a reentry inventory into an alumina emissions CSV for csv_to_alumina_events.py.

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python (numpy).

Steps:
  1. Read the ESA aluminium ablation profiles (rocket body and satellite):
       Altitude,Al Mass Loss Fraction
     Rows may be unsorted and may repeat altitudes (digitised curves). The
     profile is treated as a mass loss fraction per km of altitude and linearly
     interpolated between points. A stable sort keeps repeated altitudes in file
     order, so vertical steps in the digitised curve are kept.
  2. If a profile does not integrate to 1 over altitude, normalise it so it can
     be applied to any reentry mass.
  3. Integrate each profile over 5 km bands, 30 to 90 km. Band A holds the mass
     lost in [A, A+5) km. The 90 km band also holds everything at or above 90 km,
     and the 30 km band also holds everything below 30 km, so mass is conserved.
  4. For each reentry in the inventory, multiply its Al2O3 mass by the banded
     profile for its category:
       P (payload) and C (component)  -> satellite profile
       S0-S4 (rocket stages)           -> rocket body profile
     Al2O3 in the inventory is in tonnes; the output is in kg.
     Reentries with no Al2O3 or no reentry_lat_deg/reentry_lon_deg are skipped
     and reported.

Output CSV (one row per reentry and altitude band):
  id,start_date,lat,long,altitude_km,emitted_mass_kg
  1,2025-01-01 02:40:00,-20.8,33.9,90,0.4012

Usage:
  reentry_csv_to_alumina_emissions.py
  reentry_csv_to_alumina_emissions.py --reentries 2025_emissions.csv -o out.csv
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

DATA_DIR = (Path.home() / "Library/Mobile Documents/com~apple~CloudDocs"
            / "MIT/CARMA/Reentry Emissions")

ALT_BOTTOM_KM = 30
ALT_TOP_KM = 90
ALT_STEP_KM = 5
BANDS_KM = np.arange(ALT_TOP_KM, ALT_BOTTOM_KM - 1, -ALT_STEP_KM)  # 90, 85, ..., 30

TONNE_TO_KG = 1.0e3

SATELLITE_CATEGORIES = {"P", "C"}
ROCKET_BODY_RE = re.compile(r"^S\d+$")

DATE_RE = re.compile(r"^(\d{4}-\d{2}-\d{2})[T ](\d{2}:\d{2}(?::\d{2})?)")


class InputError(Exception):
    pass


def read_profile(path):
    """Return (altitude_km, fraction) arrays sorted by altitude."""
    alt, frac = [], []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        if len(header) < 2:
            raise InputError(f"{path}: expected 2 columns, got header {header}")
        for lineno, row in enumerate(reader, start=2):
            if not row or not "".join(row).strip():
                continue
            try:
                a, m = float(row[0]), float(row[1])
            except (ValueError, IndexError) as e:
                raise InputError(f"{path}:{lineno}: bad row {row} ({e})")
            if not (np.isfinite(a) and np.isfinite(m)) or m < 0.0:
                raise InputError(f"{path}:{lineno}: invalid altitude/fraction {row}")
            alt.append(a)
            frac.append(m)
    order = np.argsort(alt, kind="stable")
    return np.asarray(alt)[order], np.asarray(frac)[order]


def integrate(alt, frac, lo=-np.inf, hi=np.inf):
    """Exact integral over [lo, hi] km of the piecewise-linear profile.

    Nothing is extrapolated beyond the first and last profile altitudes.
    Zero-width segments (repeated altitudes) contribute nothing.
    """
    a0, a1, m0, m1 = alt[:-1], alt[1:], frac[:-1], frac[1:]
    keep = a1 > a0
    a0, a1, m0, m1 = a0[keep], a1[keep], m0[keep], m1[keep]
    x0 = np.clip(lo, a0, a1)
    x1 = np.clip(hi, a0, a1)
    slope = (m1 - m0) / (a1 - a0)
    y0 = m0 + slope * (x0 - a0)
    y1 = m0 + slope * (x1 - a0)
    return float(np.sum(0.5 * (y0 + y1) * (x1 - x0)))


def band_profile(name, alt, frac, tol):
    """Fraction of the ablated mass in each output band (ordered as BANDS_KM).

    Normalised so the bands sum to 1 if the profile's integral is not 1.
    """
    total = integrate(alt, frac)
    if total <= 0.0:
        raise InputError(f"{name}: mass loss profile integrates to {total}")
    if abs(total - 1.0) <= tol:
        print(f"{name}: profile integrates to {total:.6f} over altitude; "
              "no normalisation needed")
        total = 1.0
    else:
        print(f"{name}: profile integrates to {total:.6f} over altitude "
              f"(point sum {frac.sum():.6f}); normalising to 1")
    banded = []
    for b in BANDS_KM:
        lo = b if b > ALT_BOTTOM_KM else -np.inf
        hi = b + ALT_STEP_KM if b < ALT_TOP_KM else np.inf
        banded.append(integrate(alt, frac, lo, hi) / total)
    below = integrate(alt, frac, hi=ALT_BOTTOM_KM) / total
    above = integrate(alt, frac, lo=ALT_TOP_KM) / total
    return np.array(banded), below, above


def format_date(text, where):
    m = DATE_RE.match(text.strip())
    if not m:
        raise InputError(f"{where}: cannot parse date {text!r}")
    hms = m.group(2) if m.group(2).count(":") == 2 else m.group(2) + ":00"
    return f"{m.group(1)} {hms}"


def profile_for(category):
    if category in SATELLITE_CATEGORIES:
        return "satellite"
    if ROCKET_BODY_RE.match(category):
        return "rocket_body"
    return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--reentries", type=Path, default=DATA_DIR / "2025_emissions.csv")
    p.add_argument("--rocket-profile", type=Path,
                   default=DATA_DIR / "esa_rocket_body_ablation_profile.csv")
    p.add_argument("--satellite-profile", type=Path,
                   default=DATA_DIR / "esa_satellite_ablation_profile.csv")
    p.add_argument("-o", "--output", type=Path,
                   default=Path(__file__).with_name("alumina_2025_reentry_emissions.csv"))
    p.add_argument("--integral-tol", type=float, default=1e-6,
                   help="tolerance for 'profile already integrates to 1'")
    args = p.parse_args(argv)

    profiles = {}
    for name, path in (("rocket_body", args.rocket_profile),
                       ("satellite", args.satellite_profile)):
        alt, frac = read_profile(path)
        banded, below, above = band_profile(name, alt, frac, args.integral_tol)
        assert abs(banded.sum() - 1.0) < 1e-12
        profiles[name] = banded
        print(f"  {below:.4%} below {ALT_BOTTOM_KM} km (added to the {ALT_BOTTOM_KM} km band), "
              f"{above:.4%} at or above {ALT_TOP_KM} km")
        print("  band fractions: " + ", ".join(
            f"{b}:{f:.4f}" for b, f in zip(BANDS_KM, banded)))

    # Mass bookkeeping per category, in tonnes of Al2O3.
    stats = {}

    def tally(cat, key, mass_t):
        s = stats.setdefault(cat, {"written": [0, 0.0], "no_location": [0, 0.0],
                                   "zero_mass": [0, 0.0], "unknown_category": [0, 0.0]})
        s[key][0] += 1
        s[key][1] += mass_t

    n_rows = 0
    total_kg = 0.0
    events = []  # (start_date, lat, lon, mass_kg per band)
    with open(args.reentries, newline="") as fin:
        reader = csv.DictReader(fin)
        needed = ["date", "category", "Al2O3", "reentry_lat_deg", "reentry_lon_deg"]
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            raise InputError(f"{args.reentries}: missing column(s) {missing}")

        for lineno, r in enumerate(reader, start=2):
            where = f"{args.reentries}:{lineno}"
            cat = r["category"].strip()
            mass_t = float(r["Al2O3"] or 0.0)
            prof = profile_for(cat)
            if prof is None:
                tally(cat, "unknown_category", mass_t)
                continue
            if mass_t <= 0.0:
                tally(cat, "zero_mass", mass_t)
                continue
            lat_s, lon_s = r["reentry_lat_deg"].strip(), r["reentry_lon_deg"].strip()
            if not lat_s or not lon_s:
                tally(cat, "no_location", mass_t)
                continue
            lat, lon = float(lat_s), float(lon_s)
            if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 360.0):
                raise InputError(f"{where}: lat/lon {lat}, {lon} out of range")
            start = format_date(r["date"], where)

            events.append((start, lat_s, lon_s, mass_t * TONNE_TO_KG * profiles[prof]))
            tally(cat, "written", mass_t)

    # csv_to_alumina_events.py requires chronological order.
    events.sort(key=lambda e: e[0])
    with open(args.output, "w", newline="") as fout:
        writer = csv.writer(fout, lineterminator="\n")
        writer.writerow(["id", "start_date", "lat", "long", "altitude_km", "emitted_mass_kg"])
        for event_id, (start, lat_s, lon_s, mass_kg) in enumerate(events, start=1):
            for b, m in zip(BANDS_KM, mass_kg):
                writer.writerow([event_id, start, lat_s, lon_s, int(b), f"{m:.6g}"])
                n_rows += 1
            total_kg += mass_kg.sum()

    print(f"\nWrote {args.output}: {len(events)} reentries, {n_rows} rows, "
          f"{total_kg:,.1f} kg Al2O3")
    print(f"\n{'category':>8} {'profile':>11} | {'written':>16} | {'no location':>16} | "
          f"{'zero Al2O3':>8}   (count, tonnes Al2O3)")
    for cat in sorted(stats):
        s = stats[cat]
        prof = profile_for(cat) or "UNKNOWN"
        print(f"{cat:>8} {prof:>11} | {s['written'][0]:5d} {s['written'][1]:10.3f} | "
              f"{s['no_location'][0]:5d} {s['no_location'][1]:10.3f} | "
              f"{s['zero_mass'][0]:8d}")
        if s["unknown_category"][0]:
            print(f"         WARNING: {s['unknown_category'][0]} rows skipped "
                  f"({s['unknown_category'][1]:.3f} t), unknown category")
    dropped = sum(s["no_location"][1] for s in stats.values())
    kept = sum(s["written"][1] for s in stats.values())
    print(f"\nAl2O3 written: {kept:.3f} t; dropped for missing location: {dropped:.3f} t "
          f"({dropped / (kept + dropped):.1%} of non-zero Al2O3)")


if __name__ == "__main__":
    try:
        main()
    except InputError as e:
        sys.exit(f"error: {e}")
