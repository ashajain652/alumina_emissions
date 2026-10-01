#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Python reference for the alumina PSD remap onto CARMA bins (plan section 4).

Interpreter: /Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python

CARMA bins (base/setupbins.F90), k = vmrat, radii in cm:
  r_j    = rmin * k**((j-1)/3)
  rlow_j = r_j * (2/(k+1))**(1/3)
  rup_j  = rlow_j * k**(1/3)

Remap, assuming mass uniform in ln r inside each input bin j:
  W(i,j)       = overlap_lnr(CARMA bin i, input bin j) / dlnr_j
  psdfrac(a,i) = sum_j W(i,j) * mass_frac(a,j)
  clipped(a)   = 1 - sum_i psdfrac(a,i)     (mass outside rlow_1 .. rup_NBIN)
then each row is renormalized to sum to 1.

Usage:
  remap_psd_reference.py [alumina_test_psd.nc] [--rmin 5e-8] [--vmrat 2] [--nbin 36]
                         [--csv out.csv] [--at 67.5] [--clip-limit 0.01]
  remap_psd_reference.py --selftest
"""

import argparse
import sys
from pathlib import Path

import numpy as np

CM_TO_NM = 1.0e7
DEFAULT_PSD = Path(__file__).resolve().parent / "alumina_test_psd.nc"


def carma_bin_edges(rmin_cm=5e-8, vmrat=2.0, nbin=36):
    """Return (r, rlow, rup) in nm, following base/setupbins.F90."""
    k = float(vmrat)
    j = np.arange(1, nbin + 1)
    r = rmin_cm * k ** ((j - 1) / 3.0)
    rlow = r * (2.0 / (k + 1.0)) ** (1.0 / 3.0)
    rup = rlow * k ** (1.0 / 3.0)
    return r * CM_TO_NM, rlow * CM_TO_NM, rup * CM_TO_NM


def overlap_weights(r_edge_in, rlow, rup):
    """W[i, j]: fraction of input bin j's ln r width that lies in CARMA bin i."""
    lin = np.log(np.asarray(r_edge_in, dtype=float))
    lo_in, hi_in = lin[:-1], lin[1:]
    lo_c, hi_c = np.log(rlow), np.log(rup)
    ov = (np.minimum(hi_c[:, None], hi_in[None, :])
          - np.maximum(lo_c[:, None], lo_in[None, :]))
    ov = np.clip(ov, 0.0, None)
    return ov / (hi_in - lo_in)[None, :]


def remap(mass_frac, r_edge_in, rlow, rup):
    """Return (psdfrac renormalized [nalt, nbin], clipped fraction [nalt])."""
    mass_frac = np.atleast_2d(np.asarray(mass_frac, dtype=float))
    W = overlap_weights(r_edge_in, rlow, rup)
    raw = mass_frac @ W.T                      # [nalt, NBIN]
    total_in = mass_frac.sum(axis=1)
    kept = raw.sum(axis=1)
    clipped = 1.0 - kept / total_in
    psdfrac = raw / kept[:, None]
    return psdfrac, clipped


def interp_binfrac(psd_alt, psdfrac, z_km):
    """Linear interpolation of psdfrac rows in altitude, clamped outside the range."""
    psd_alt = np.asarray(psd_alt, dtype=float)
    if psd_alt.size == 1 or z_km <= psd_alt[0]:
        return psdfrac[0].copy()
    if z_km >= psd_alt[-1]:
        return psdfrac[-1].copy()
    a = int(np.searchsorted(psd_alt, z_km, side="right")) - 1
    w = (z_km - psd_alt[a]) / (psd_alt[a + 1] - psd_alt[a])
    return (1.0 - w) * psdfrac[a] + w * psdfrac[a + 1]


def read_psd_nc(path):
    from netCDF4 import Dataset
    with Dataset(path) as nc:
        alt = np.array(nc["psd_alt"][:], dtype=float)
        r_edge = np.array(nc["r_edge"][:], dtype=float)
        mf = np.array(nc["mass_frac"][:], dtype=float)
    if np.any(np.diff(alt) <= 0):
        raise ValueError("psd_alt must be strictly increasing")
    if np.any(np.diff(r_edge) <= 0):
        raise ValueError("r_edge must be strictly increasing")
    return alt, r_edge, mf


def print_table(alt, psdfrac, clipped, r, rlow, rup, tol=0.0):
    nz = np.where(np.any(psdfrac > tol, axis=0))[0]
    print("\nCARMA bins receiving mass (radius nm):")
    for i in nz:
        print(f"  bin {i + 1:2d}: rlow {rlow[i]:9.3f}  r {r[i]:9.3f}  rup {rup[i]:9.3f}")
    head = "alt_km  clipped   " + " ".join(f"   bin{i + 1:02d}  " for i in nz)
    print("\npsdfrac (renormalized; nonzero CARMA bins only)")
    print(head)
    for a, z in enumerate(alt):
        vals = " ".join(f"{psdfrac[a, i]:10.3e}" for i in nz)
        print(f"{z:6.2f}  {clipped[a]:8.2e}  {vals}")
    return nz


def write_csv(path, alt, psdfrac, clipped):
    nbin = psdfrac.shape[1]
    with open(path, "w") as f:
        f.write("altitude_km,clipped," + ",".join(f"bin{i + 1:02d}" for i in range(nbin)) + "\n")
        for a, z in enumerate(alt):
            f.write(f"{z:g},{clipped[a]:.10e}," +
                    ",".join(f"{v:.10e}" for v in psdfrac[a]) + "\n")


def selftest():
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("PASS  " if cond else "FAIL  ") + msg)
        ok = ok and bool(cond)

    for rmin, vmrat, nbin in ((5e-8, 2.0, 36), (2e-8, 2.0, 28), (1e-7, 3.7, 20)):
        r, rlow, rup = carma_bin_edges(rmin, vmrat, nbin)
        tag = f"rmin={rmin:g} vmrat={vmrat:g} nbin={nbin}"
        gap = np.max(np.abs(rup[:-1] - rlow[1:]) / rlow[1:])
        check(gap < 1e-12, f"CARMA edges contiguous ({tag}), max rel gap {gap:.1e}")
        check(np.all((rlow < r) & (r < rup)), f"r_j inside [rlow_j, rup_j] ({tag})")

        # Identity: PSD built on CARMA's own bins (a sub-range) maps onto itself.
        rng = np.random.default_rng(1)
        j0, j1 = 3, nbin - 4
        edges = np.concatenate([rlow[j0:j1], [rup[j1 - 1]]])
        mf = rng.random((4, j1 - j0))
        mf /= mf.sum(axis=1, keepdims=True)
        pf, cl = remap(mf, edges, rlow, rup)
        expect = np.zeros((4, nbin))
        expect[:, j0:j1] = mf
        err = np.max(np.abs(pf - expect))
        check(err < 1e-12, f"own-bin PSD maps to identity ({tag}), max err {err:.1e}")
        check(np.max(np.abs(cl)) < 1e-12, f"no clipping for own-bin PSD ({tag})")
        check(np.max(np.abs(pf.sum(axis=1) - 1)) < 1e-12, f"rows sum to 1 ({tag})")

    # Arbitrary PSD, partly outside the CARMA range: clipping and renormalization.
    r, rlow, rup = carma_bin_edges()
    edges = np.geomspace(rlow[0] * 0.5, rup[-1] * 2.0, 50)
    mf = np.ones((2, 49)) / 49.0
    pf, cl = remap(mf, edges, rlow, rup)
    check(np.max(np.abs(pf.sum(axis=1) - 1)) < 1e-12, "rows sum to 1 after clipping")
    lnw = np.log(edges[-1] / edges[0])
    exp_clip = (np.log(rlow[0] / edges[0]) + np.log(edges[-1] / rup[-1])) / lnw
    check(np.allclose(cl, exp_clip, rtol=1e-10),
          f"clipped fraction {cl[0]:.6f} equals analytic {exp_clip:.6f}")

    # Interpolation: exact at nodes, linear between, clamped outside.
    alt = np.array([45.0, 65.0, 85.0])
    tab = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]])
    check(np.allclose(interp_binfrac(alt, tab, 65.0), [0, 1]), "interp exact at node")
    check(np.allclose(interp_binfrac(alt, tab, 55.0), [0.5, 0.5]), "interp linear midpoint")
    check(np.allclose(interp_binfrac(alt, tab, 10.0), [1, 0]), "interp clamps below")
    check(np.allclose(interp_binfrac(alt, tab, 120.0), [0.5, 0.5]), "interp clamps above")

    print("\nSELFTEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("psd", nargs="?", default=str(DEFAULT_PSD), help="PSD netCDF file")
    p.add_argument("--rmin", type=float, default=5e-8, help="CARMA rmin in cm (default 5e-8)")
    p.add_argument("--vmrat", type=float, default=2.0, help="CARMA vmrat (default 2)")
    p.add_argument("--nbin", type=int, default=36, help="CARMA NBIN (default 36)")
    p.add_argument("--csv", help="write the full psdfrac table to this CSV")
    p.add_argument("--at", type=float, action="append", default=[],
                   help="also print the interpolated binfrac at this altitude (km); repeatable")
    p.add_argument("--clip-limit", type=float, default=0.01,
                   help="clipped fraction that the Fortran treats as fatal (default 0.01)")
    p.add_argument("--selftest", action="store_true", help="run the self-test and exit")
    args = p.parse_args(argv)

    if args.selftest:
        return selftest()

    r, rlow, rup = carma_bin_edges(args.rmin, args.vmrat, args.nbin)
    alt, r_edge, mf = read_psd_nc(args.psd)
    psdfrac, clipped = remap(mf, r_edge, rlow, rup)

    print(f"PSD file: {args.psd}")
    print(f"CARMA: rmin={args.rmin:g} cm, vmrat={args.vmrat:g}, NBIN={args.nbin}; "
          f"radius range {rlow[0]:.4g} .. {rup[-1]:.4g} nm")
    print(f"PSD:   {len(alt)} altitudes, {len(r_edge) - 1} bins, "
          f"radius {r_edge[0]:.4g} .. {r_edge[-1]:.4g} nm")
    print("\nClipped fraction per altitude:")
    for z, c in zip(alt, clipped):
        flag = "  ** exceeds limit (endrun unless carma_emis_psd_allow_clip)" \
            if c > args.clip_limit else ""
        print(f"  {z:7.2f} km  {c:.3e}{flag}")
    print_table(alt, psdfrac, clipped, r, rlow, rup)

    for z in args.at:
        bf = interp_binfrac(alt, psdfrac, z)
        nz = np.where(bf > 0)[0]
        print(f"\nbinfrac at {z:g} km: " +
              ", ".join(f"bin{i + 1}={bf[i]:.4e}" for i in nz) + f"  (sum {bf.sum():.12f})")

    if args.csv:
        write_csv(args.csv, alt, psdfrac, clipped)
        print(f"\nWrote {args.csv}")

    return 1 if np.any(clipped > args.clip_limit) else 0


if __name__ == "__main__":
    sys.exit(main())
