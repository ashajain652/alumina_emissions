#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Convert an alumina PSD CSV into the CARMA alumina PSD netCDF file.

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
(numpy + netCDF4; the miniforge base env has neither).

Input CSV (see ALUMINA_CARMA_PLAN.md section 4):
  altitude_km,bin,d_low_um,d_center_um,d_high_um,dNdlnr_cm-3
Diameters in um; dNdlnr is the bin-averaged dN/dlnr (= dN/dlnD).

Output netCDF (classic 64-bit offset format):
  dims psd_alt, psd_bin, psd_bin_edge (= psd_bin + 1)
  psd_alt(psd_alt)            km above ground, increasing
  r_edge(psd_bin_edge)        radius in nm (= diameter um * 500), increasing
  mass_frac(psd_alt, psd_bin) each row sums to 1

Mass per bin, with dN/dlnr constant in ln r inside the bin:
  m_i ∝ dNdlnr_i * (r_hi^3 - r_lo^3) / 3, normalized per altitude.
Density cancels in the normalized fractions.

Usage:
  csv_to_alumina_psd.py alumina_test_psd.csv alumina_test_psd.nc
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

UM_DIAM_TO_NM_RADIUS = 500.0
EDGE_RTOL = 1e-6  # CSV edges are printed with 6 decimals

REQUIRED_COLUMNS = ["altitude_km", "bin", "d_low_um", "d_high_um", "dNdlnr_cm-3"]


class InputError(Exception):
    pass


def close(a, b, rtol=EDGE_RTOL):
    return abs(a - b) <= rtol * max(abs(a), abs(b))


def read_psd_csv(path):
    by_alt = defaultdict(dict)
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        header = [h.strip() for h in (reader.fieldnames or [])]
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise InputError(f"{path}: missing column(s) {missing}; found {header}")
        for lineno, raw in enumerate(reader, start=2):
            r = {k.strip(): (v.strip() if v is not None else "") for k, v in raw.items()}
            if not any(r.values()):
                continue
            where = f"{path}:{lineno}"
            try:
                alt = float(r["altitude_km"])
                b = int(r["bin"])
                dlo = float(r["d_low_um"])
                dhi = float(r["d_high_um"])
                val = float(r["dNdlnr_cm-3"])
            except ValueError as e:
                raise InputError(f"{where}: bad number ({e})")
            if not all(math.isfinite(x) for x in (alt, dlo, dhi, val)):
                raise InputError(f"{where}: non-finite value")
            if not 0.0 < dlo < dhi:
                raise InputError(f"{where}: need 0 < d_low < d_high (got {dlo}, {dhi})")
            if val < 0.0:
                raise InputError(f"{where}: negative dNdlnr {val}")
            if b in by_alt[alt]:
                raise InputError(f"{where}: duplicate bin {b} at altitude {alt} km")
            by_alt[alt][b] = (dlo, dhi, val)
    if not by_alt:
        raise InputError(f"{path}: no data rows")
    return by_alt


def build(by_alt, path):
    alts = sorted(by_alt)  # increasing altitude
    ref_alt = alts[0]
    ref_bins = sorted(by_alt[ref_alt])
    ref = [by_alt[ref_alt][b][:2] for b in ref_bins]  # sorted by bin number

    # Order bins by lower edge, check contiguity.
    order = sorted(range(len(ref)), key=lambda i: ref[i][0])
    ref_sorted = [ref[i] for i in order]
    bins_sorted = [ref_bins[i] for i in order]
    for (lo0, hi0), (lo1, hi1), b1 in zip(ref_sorted, ref_sorted[1:], bins_sorted[1:]):
        if not close(hi0, lo1):
            raise InputError(f"{path}: bins not contiguous: d_high {hi0} != next "
                             f"d_low {lo1} (bin {b1})")

    d_edges = [lo for lo, _ in ref_sorted] + [ref_sorted[-1][1]]

    nalt, nbin = len(alts), len(bins_sorted)
    dndlnr = np.zeros((nalt, nbin))
    for a, alt in enumerate(alts):
        bins = by_alt[alt]
        if sorted(bins) != sorted(bins_sorted):
            raise InputError(f"{path}: altitude {alt} km has bins {sorted(bins)}, "
                             f"altitude {ref_alt} km has {sorted(bins_sorted)}")
        for i, b in enumerate(bins_sorted):
            dlo, dhi, val = bins[b]
            rlo, rhi = ref_sorted[i]
            if not (close(dlo, rlo) and close(dhi, rhi)):
                raise InputError(f"{path}: bin {b} edges at {alt} km ({dlo}, {dhi}) "
                                 f"differ from those at {ref_alt} km ({rlo}, {rhi})")
            dndlnr[a, i] = val

    r_edge = np.array(d_edges) * UM_DIAM_TO_NM_RADIUS  # nm radius
    # Force exact contiguity using the shared edge array.
    vol = (r_edge[1:] ** 3 - r_edge[:-1] ** 3) / 3.0
    mass = dndlnr * vol[None, :]
    tot = mass.sum(axis=1)
    if np.any(tot <= 0.0):
        bad = [alts[a] for a in np.where(tot <= 0.0)[0]]
        raise InputError(f"{path}: zero total mass at altitude(s) {bad} km")
    mass_frac = mass / tot[:, None]
    return np.array(alts), r_edge, mass_frac


def write_nc(out, alts, r_edge, mass_frac, source):
    with Dataset(out, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.createDimension("psd_alt", len(alts))
        nc.createDimension("psd_bin", len(r_edge) - 1)
        nc.createDimension("psd_bin_edge", len(r_edge))
        v = nc.createVariable("psd_alt", "f8", ("psd_alt",))
        v.units = "km"
        v.long_name = "altitude above ground of PSD rows"
        v[:] = alts
        v = nc.createVariable("r_edge", "f8", ("psd_bin_edge",))
        v.units = "nm"
        v.long_name = "PSD bin edge radius"
        v[:] = r_edge
        v = nc.createVariable("mass_frac", "f8", ("psd_alt", "psd_bin"))
        v.units = "1"
        v.long_name = "fraction of emitted mass in each PSD bin (rows sum to 1)"
        v[:] = mass_frac
        nc.title = "CARMA alumina emission particle size distribution"
        nc.source_csv = str(source)
        nc.history = "created by csv_to_alumina_psd.py"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="input PSD CSV")
    p.add_argument("out", help="output PSD netCDF")
    args = p.parse_args(argv)

    src = Path(args.csv).resolve()
    try:
        alts, r_edge, mass_frac = build(read_psd_csv(src), src)
    except InputError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    write_nc(args.out, alts, r_edge, mass_frac, src)

    print(f"Wrote {args.out}: {len(alts)} altitudes x {len(r_edge) - 1} bins")
    print(f"  psd_alt = {', '.join(f'{a:g}' for a in alts)} km")
    print(f"  r_edge  = {r_edge[0]:.4g} .. {r_edge[-1]:.4g} nm")
    print(f"  max |row sum - 1| = {np.abs(mass_frac.sum(axis=1) - 1).max():.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
