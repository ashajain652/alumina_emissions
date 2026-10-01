#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Build the per-band Al2O3 refractive index table for the CARMA alumina model.

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python

Inputs (two columns, no header, wavelength in um; points may be unsorted,
since both files were digitized from published plots):
  tropf_real_alumina_2.csv                real part n
  tropf_kapoor_combined_img_alumina.csv   imaginary part k

CARMA does one Mie calculation per RRTMG band, at the band-centre wavelength,
with one refractive index per band. The band centres are computed here exactly
as in carma_intr.F90 (CARMA_Register): the arithmetic midpoint in wavelength of
the band edges from radconstants.F90. Order is the CARMA wave order: the 16 LW
bands, then the 14 SW bands.

n and k are sampled at the band centre (no band averaging):
  n: linear in ln(lambda)
  k: linear in ln(lambda) and ln(k)
The data cover 0.0096-642 um, so every band centre is inside the data range and
no extrapolation is done (the script stops if one is not).

Outputs:
  alumina_refidx_bands.csv   band table (edges, centre, n, k)
  stdout                     Fortran table to paste into
                             models/alumina/carma_model_mod.F90

Usage:
  make_alumina_refidx.py [--real FILE] [--imag FILE] [--out CSV]
"""

import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# RRTMG band limits (cm-1), from src/physics/rrtmg/radconstants.F90.
LW_LO = [10, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2390, 2600]
LW_HI = [350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2390, 2600, 3250]
SW_LO = [2600, 3250, 4000, 4650, 5150, 6150, 7700, 8050, 12850, 16000, 22650, 29000, 38000, 820]
SW_HI = [3250, 4000, 4650, 5150, 6150, 7700, 8050, 12850, 16000, 22650, 29000, 38000, 50000, 2600]


def load(path):
    a = np.loadtxt(path, delimiter=",")
    a = a[np.argsort(a[:, 0])]
    if np.any(np.diff(a[:, 0]) <= 0):
        raise SystemExit(f"{path}: repeated wavelengths")
    if np.any(a[:, 1] <= 0):
        raise SystemExit(f"{path}: non-positive index values")
    return a


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--real", default=HERE / "tropf_real_alumina_2.csv")
    p.add_argument("--imag", default=HERE / "tropf_kapoor_combined_img_alumina.csv")
    p.add_argument("--out", default=HERE / "alumina_refidx_bands.csv")
    args = p.parse_args()

    re = load(args.real)
    im = load(args.imag)

    wn_lo = np.array(LW_LO + SW_LO, float)
    wn_hi = np.array(LW_HI + SW_HI, float)
    lam_lo = 1e4 / wn_hi  # um
    lam_hi = 1e4 / wn_lo
    cen = lam_lo + (lam_hi - lam_lo) / 2.0  # as carma_intr.F90
    names = [f"LW{i+1:02d}" for i in range(16)] + [f"SW{i+1:02d}" for i in range(14)]

    for a, f in ((re, args.real), (im, args.imag)):
        if cen.min() < a[0, 0] or cen.max() > a[-1, 0]:
            raise SystemExit(f"{f}: band centres {cen.min():.4g}-{cen.max():.4g} um outside data range")

    n = np.interp(np.log(cen), np.log(re[:, 0]), re[:, 1])
    k = np.exp(np.interp(np.log(cen), np.log(im[:, 0]), np.log(im[:, 1])))

    with open(args.out, "w") as fh:
        fh.write("band,wave_lo_um,wave_hi_um,wave_center_um,n,k\n")
        for i in range(30):
            fh.write(f"{names[i]},{lam_lo[i]:.6f},{lam_hi[i]:.6f},{cen[i]:.6f},{n[i]:.6f},{k[i]:.6e}\n")
    print(f"! wrote {args.out}")
    print(f"! source: {Path(args.real).name} (n), {Path(args.imag).name} (k)")

    # Fortran table
    print("    refidx(:) = (/ &")
    for i in range(30):
        sep = ", &" if i < 29 else " /)"
        print(f"      ({n[i]:.6f}_f, {k[i]:.6e}_f){sep}   ! {names[i]} {cen[i]:9.4f} um")


if __name__ == "__main__":
    main()
