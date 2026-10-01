#!/usr/bin/env python3
"""Generate test inputs for the CARMA alumina model.

Writes, next to this script:
  alumina_test_emissions.csv  10 reentry emission events (dirac pulses) over 2 days
  alumina_test_psd.csv        lognormal dN/dlnr on diameter bins, 45-85 km

Only the standard library is used. The random seed is fixed, so rerunning
reproduces the files exactly.
"""

import csv
import math
import random
from datetime import datetime, timedelta
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent
SEED = 42

# Altitude grid shared by both files (km above ground, not sea level).
ALTITUDES_KM = list(range(85, 40, -5))  # 85, 80, ..., 45

# ---------------------------------------------------------------------------
# Emission events
# ---------------------------------------------------------------------------
# The case FWsc2010climo.CARMA_dev starts at 0001-01-01 00:00 (noleap calendar)
# and runs for 2 days. The events are spread evenly over that window, each at the
# centre of its slot: (i + 0.5) * 48 h / 10.
RUN_START = datetime(1, 1, 1, 0, 0, 0)
RUN_HOURS = 48.0
N_EVENTS = 10

# Vertical shape of the ablation mass: a Gaussian in altitude, normalized over
# the 9 levels. Test values only.
PROFILE_PEAK_KM = 70.0
PROFILE_SIGMA_KM = 8.0

# Total alumina mass per event (kg), drawn uniformly from this range.
EVENT_MASS_RANGE_KG = (50.0, 300.0)


def make_emissions(rng):
    weights = [math.exp(-0.5 * ((z - PROFILE_PEAK_KM) / PROFILE_SIGMA_KM) ** 2)
               for z in ALTITUDES_KM]
    wsum = sum(weights)
    weights = [w / wsum for w in weights]

    rows = []
    slot_hours = RUN_HOURS / N_EVENTS
    for i in range(N_EVENTS):
        start = RUN_START + timedelta(hours=(i + 0.5) * slot_hours)
        # Uniform on the sphere: lat = asin(u), u uniform in [-1, 1].
        lat = math.degrees(math.asin(rng.uniform(-1.0, 1.0)))
        lon = rng.uniform(-180.0, 180.0)
        total = rng.uniform(*EVENT_MASS_RANGE_KG)
        for z, w in zip(ALTITUDES_KM, weights):
            rows.append({
                "id": i + 1,
                # Format by hand: strftime's %Y does not zero-pad year 1 on every platform.
                "start_date": (f"{start.year:04d}-{start.month:02d}-{start.day:02d} "
                               f"{start.hour:02d}:{start.minute:02d}:{start.second:02d}"),
                "lat": f"{lat:.4f}",
                "long": f"{lon:.4f}",
                "altitude_km": z,
                "emitted_mass_kg": f"{total * w:.4f}",
            })

    # Chronological, latest at the bottom; within an event, altitude 85 -> 45 km.
    rows.sort(key=lambda r: (r["start_date"], -r["altitude_km"]))

    path = OUT_DIR / "alumina_test_emissions.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path, rows


# ---------------------------------------------------------------------------
# Particle size distribution
# ---------------------------------------------------------------------------
# Number lognormal in DIAMETER, truncated to 0.1-1 um. The median is the
# geometric midpoint of the range. dN/dlnr = dN/dlnD, because ln r = ln D - ln 2.
D_MIN_UM, D_MAX_UM = 0.1, 1.0
N_BINS = 20
SIGMA_G = 1.6
D_MEDIAN_UM = math.sqrt(D_MIN_UM * D_MAX_UM)
N_TOTAL_CM3 = 1.0  # Arbitrary: the model only uses the shape.


def lognormal_cdf(d):
    return 0.5 * (1.0 + math.erf(math.log(d / D_MEDIAN_UM) / (math.sqrt(2.0) * math.log(SIGMA_G))))


def make_psd():
    ratio = (D_MAX_UM / D_MIN_UM) ** (1.0 / N_BINS)
    edges = [D_MIN_UM * ratio ** k for k in range(N_BINS + 1)]
    edges[-1] = D_MAX_UM  # Avoid float drift at the top edge.

    bins = []
    for k in range(N_BINS):
        lo, hi = edges[k], edges[k + 1]
        n_bin = N_TOTAL_CM3 * (lognormal_cdf(hi) - lognormal_cdf(lo))
        # Bin-averaged dN/dlnD, so that sum(dNdlnr * dlnD) is exactly the number in range.
        bins.append((k + 1, lo, math.sqrt(lo * hi), hi, n_bin / math.log(hi / lo)))

    rows = []
    for z in ALTITUDES_KM:
        for b, lo, ctr, hi, val in bins:
            rows.append({
                "altitude_km": z,
                "bin": b,
                "d_low_um": f"{lo:.6f}",
                "d_center_um": f"{ctr:.6f}",
                "d_high_um": f"{hi:.6f}",
                "dNdlnr_cm-3": f"{val:.6e}",
            })

    path = OUT_DIR / "alumina_test_psd.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    in_range = lognormal_cdf(D_MAX_UM) - lognormal_cdf(D_MIN_UM)
    return path, in_range


def main():
    rng = random.Random(SEED)
    epath, erows = make_emissions(rng)
    ppath, in_range = make_psd()

    total_kg = sum(float(r["emitted_mass_kg"]) for r in erows)
    print(f"Wrote {epath} ({len(erows)} rows, {N_EVENTS} events)")
    print(f"  total emitted mass = {total_kg:.4f} kg")
    print(f"Wrote {ppath} ({len(ALTITUDES_KM)} altitudes x {N_BINS} bins)")
    print(f"  fraction of the untruncated lognormal (by number) inside "
          f"{D_MIN_UM}-{D_MAX_UM} um = {in_range:.4f}")


if __name__ == "__main__":
    main()
