#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Write a CARMA alumina PSD netCDF file for a single particle size at all altitudes.

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python

The file has the layout that alumina_psd_mod.F90 reads (plan section 4):
  dims psd_alt = 1, psd_bin = 1, psd_bin_edge = 2
  psd_alt(psd_alt)            km above ground (a single row, used at every altitude)
  r_edge(psd_bin_edge)        radius in nm, a very narrow bin around the given size
  mass_frac(psd_alt, psd_bin) = 1

alumina_psd_binfrac clamps outside of the psd_alt range, so a single row applies
to every emission altitude. The remap splits the input bin among CARMA bins by
ln r overlap; since the input bin is very narrow, all of the mass goes into the
one CARMA bin that contains the given size. The emitted particles therefore have
the size of that CARMA bin, not exactly the size given here.

The size is a diameter in um by default, as in the PSD CSV; use --radius to give
a radius instead.

--no-check skips the check against the CARMA bins, so a size outside the default
grid can be written for a model built with a wider grid.

Usage:
  make_monodisperse_psd.py 0.5                       # 0.5 um diameter
  make_monodisperse_psd.py 0.25 --radius -o alumina_r0.25um_psd.nc
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from remap_psd_reference import carma_bin_edges, remap

UM_TO_NM = 1000.0
HALF_WIDTH = 1.0e-6  # half width of the input bin in ln r


def write_nc(out, r_nm, r_edge, size_um, is_radius):
    with Dataset(out, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.createDimension("psd_alt", 1)
        nc.createDimension("psd_bin", 1)
        nc.createDimension("psd_bin_edge", 2)
        v = nc.createVariable("psd_alt", "f8", ("psd_alt",))
        v.units = "km"
        v.long_name = "altitude above ground of PSD rows"
        v[:] = [0.0]
        v = nc.createVariable("r_edge", "f8", ("psd_bin_edge",))
        v.units = "nm"
        v.long_name = "PSD bin edge radius"
        v[:] = r_edge
        v = nc.createVariable("mass_frac", "f8", ("psd_alt", "psd_bin"))
        v.units = "1"
        v.long_name = "fraction of emitted mass in each PSD bin (rows sum to 1)"
        v[:] = [[1.0]]
        nc.title = "CARMA alumina emission particle size distribution (monodisperse)"
        nc.particle_radius_nm = r_nm
        nc.input_size = f"{size_um:g} um {'radius' if is_radius else 'diameter'}"
        nc.history = "created by make_monodisperse_psd.py"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("size_um", type=float, help="particle size in um (diameter unless --radius)")
    p.add_argument("--radius", action="store_true", help="size_um is a radius, not a diameter")
    p.add_argument("-o", "--out", help="output netCDF (default: alumina_mono_d<size>um_psd.nc)")
    p.add_argument("--rmin", type=float, default=5e-8, help="CARMA rmin (cm), default 5e-8")
    p.add_argument("--vmrat", type=float, default=2.0, help="CARMA vmrat, default 2")
    p.add_argument("--nbin", type=int, default=36, help="CARMA NBIN, default 36")
    p.add_argument("--no-check", action="store_true",
                   help="skip the check against the CARMA bins (rmin, vmrat, nbin)")
    args = p.parse_args(argv)

    if not args.size_um > 0.0:
        print(f"ERROR: size must be > 0 (got {args.size_um})", file=sys.stderr)
        return 1

    r_nm = args.size_um * UM_TO_NM * (1.0 if args.radius else 0.5)
    r_edge = r_nm * np.exp([-HALF_WIDTH, HALF_WIDTH])

    # Check where the size lands in the CARMA bins, with the same remap as the model.
    if not args.no_check:
        r_c, rlow, rup = carma_bin_edges(args.rmin, args.vmrat, args.nbin)
        if r_nm <= rlow[0] or r_nm >= rup[-1]:
            print(f"ERROR: radius {r_nm:g} nm is outside the CARMA bins "
                  f"({rlow[0]:.4g} .. {rup[-1]:.4g} nm radius)", file=sys.stderr)
            return 1
        psdfrac, clipped = remap([[1.0]], r_edge, rlow, rup)
        bins = np.nonzero(psdfrac[0] > 0.0)[0]
        if len(bins) > 1:
            print(f"WARNING: radius {r_nm:g} nm is on a CARMA bin edge; the mass is split "
                  f"between bins {', '.join(str(b + 1) for b in bins)}", file=sys.stderr)

    if args.out:
        out = Path(args.out)
    else:
        tag = f"{'r' if args.radius else 'd'}{args.size_um:g}um"
        out = Path(__file__).resolve().parent / f"alumina_mono_{tag}_psd.nc"
    write_nc(out, r_nm, r_edge, args.size_um, args.radius)

    print(f"Wrote {out}")
    print(f"  input size   = {args.size_um:g} um {'radius' if args.radius else 'diameter'}"
          f" -> radius {r_nm:.6g} nm, used at all altitudes")
    if args.no_check:
        print("  CARMA bin check skipped (--no-check)")
        return 0
    for b in bins:
        print(f"  CARMA bin {b + 1:2d}: rlow {rlow[b]:.4g} nm, r {r_c[b]:.4g} nm, "
              f"rup {rup[b]:.4g} nm, mass fraction {psdfrac[0, b]:.6f}")
    print(f"  clipped fraction = {clipped[0]:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
