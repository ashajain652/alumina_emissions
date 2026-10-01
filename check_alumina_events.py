#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Print the event list and total mass from a CARMA alumina event netCDF file.

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python

Usage:
  check_alumina_events.py [alumina_test_events.nc] [--expect 1693.9828]
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

DEFAULT = Path(__file__).resolve().parent / "alumina_test_events.nc"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("events", nargs="?", default=str(DEFAULT), help="event netCDF file")
    p.add_argument("--expect", type=float,
                   help="expected total mass in kg; exit 1 if it differs by > 1e-9 relative")
    args = p.parse_args(argv)

    with Dataset(args.events) as nc:
        eid = np.array(nc["event_id"][:])
        date = np.array(nc["date"][:])
        sec = np.array(nc["datesec"][:])
        lat = np.array(nc["lat"][:], dtype=float)
        lon = np.array(nc["lon"][:], dtype=float)
        alt = np.array(nc["altitude"][:], dtype=float)
        mass = np.array(nc["mass"][:], dtype=float)
        attrs = {k: nc.getncattr(k) for k in nc.ncattrs()}

    print(f"Event file: {args.events}")
    for k in ("source_csv", "calendar", "total_mass_kg"):
        print(f"  {k} = {attrs.get(k, '(missing)')}")
    print(f"  rows = {len(eid)}")

    print(f"\n{'id':>4} {'date':>8} {'time':>8} {'lat':>9} {'lon':>9} "
          f"{'nlev':>4} {'alt_km':>11} {'mass_kg':>12}")
    order = []
    for e in eid:
        if e not in order:
            order.append(e)
    for e in order:
        m = eid == e
        i = np.where(m)[0][0]
        s = int(sec[i])
        hhmmss = f"{s // 3600:02d}:{s % 3600 // 60:02d}:{s % 60:02d}"
        print(f"{e:4d} {date[i]:08d} {hhmmss:>8} {lat[i]:9.4f} {lon[i]:9.4f} "
              f"{m.sum():4d} {alt[m].min():5.1f}-{alt[m].max():<5.1f} "
              f"{math.fsum(mass[m]):12.4f}")

    total = math.fsum(mass)
    print(f"\n{len(order)} events, total mass = {total:.6f} kg")

    rc = 0
    attr_total = attrs.get("total_mass_kg")
    if attr_total is not None and abs(float(attr_total) - total) > 1e-9 * max(total, 1.0):
        print(f"WARNING: total_mass_kg attribute {attr_total} != sum of mass {total}")
        rc = 1
    if args.expect is not None:
        rel = abs(total - args.expect) / max(abs(args.expect), 1e-30)
        status = "OK" if rel <= 1e-9 else "MISMATCH"
        print(f"expected {args.expect} kg: relative difference {rel:.2e} -> {status}")
        if status != "OK":
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
