#!/usr/bin/env python3
"""Write one CARMA alumina PSD netCDF file per droplet size from the droplet final radii.

Interpreter: any python3 with numpy + scipy (written with scipy.io.netcdf_file).

Inputs, from run_droplet_emission_profiles.py (Droplet Emission Profiles/):
  - droplet_<D>um_final_radius_profile.csv:
      Altitude,Surviving Al Mass Fraction,r_mass_mean_um,r_min_um,r_max_um
    Per 1 km final-altitude bin (the same bins and altitudes as the emission
    files): the mass-averaged final Al droplet radius (sum m r / sum m over the
    release parcels landing there) and the min and max radius; empty where no
    mass lands.
  - droplet_emission_summary.csv: rho_droplet_kg_m3 for each size.

Steps for each droplet size:
  1. Spread test, per altitude: if ln(r_max / r_min) in a bin is narrower than
     one CARMA bin, ln(vmrat) / 3, that altitude is treated as a single size.
     Otherwise the script stops: a polydisperse row would be needed.
  2. The droplets cool and oxidise to alumina. Each Al droplet of mass m
     becomes Al2O3 of mass m * M_Al2O3 / (2 M_Al), so the radius scales by
       (M_Al2O3 / (2 M_Al) * rho_droplet / rho_alumina)^(1/3).
  3. Write a PSD with one row per altitude that has mass (psd_alt = the bin
     centres, so every emission row falls exactly on a PSD row and the model's
     linear interpolation in altitude never blends two rows). Each row is
     monodisperse at that altitude's mass-averaged alumina radius: r_edge holds
     a very narrow bin around each distinct radius, with zero-mass gap bins in
     between, and the row's mass_frac is 1 in its own narrow bin.

Placement in the CARMA bins is reported for the grid of the matching
alumina_d<D>um model (rmin = d0/2, vmrat 2, 16 bins; the 25, 35 and 45 um models
do not exist yet), for the planned grid with rmin one bin lower
(d0/2 * 2^(-1/3); see SIZE_MODEL_RMIN_PLAN.md) and, if it fits, the base alumina
grid (rmin 5e-8 cm, 36 bins),
as the share of the emitted mass whose radius falls in each CARMA bin or outside
the grid (clipped; the model stops if a row is more than 1 % clipped).

Output: <outdir>/alumina_droplet_d<D>um_psd.nc

Usage:
  droplet_summary_to_alumina_psd.py
  droplet_summary_to_alumina_psd.py --rho-alumina 3.6 -o other_dir
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

from droplet_profiles_to_alumina_events import AL_PER_AL2O3, size_tag
from remap_psd_reference import carma_bin_edges

HERE = Path(__file__).resolve().parent
PROFILE_DIR = HERE.parent / "Droplet Emission Profiles"
RADIUS_RE = re.compile(r"^droplet_([0-9.]+)um_final_radius_profile\.csv$")

UM_TO_NM = 1000.0
HALF_WIDTH = 1.0e-6  # half width of each narrow input bin in ln r, as in make_monodisperse_psd.py
MERGE_LNR = 4.0 * HALF_WIDTH  # radii closer than this share one narrow bin
RADIUS_RTOL = 1e-12  # mass mean vs min/max of a bin, rounding

SIZE_MODEL_NBIN = 16
PLANNED_RMIN_FACTOR = 2.0 ** (-1.0 / 3.0)  # size models with rmin one bin below r0 (SIZE_MODEL_RMIN_PLAN.md)
BASE_RMIN_CM, BASE_NBIN = 5e-8, 36


class InputError(Exception):
    pass


def read_radius_profile(path):
    """Rows with mass: (altitude_km, mass_frac_per_km, r_mean_um, r_min_um, r_max_um)."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        needed = ["Altitude", "Surviving Al Mass Fraction", "r_mass_mean_um", "r_min_um",
                  "r_max_um"]
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            raise InputError(f"{path}: missing column(s) {missing}")
        for lineno, r in enumerate(reader, start=2):
            frac = float(r["Surviving Al Mass Fraction"])
            if frac <= 0.0:
                continue
            vals = [float(r[k]) for k in ("r_mass_mean_um", "r_min_um", "r_max_um")]
            if not all(math.isfinite(v) and v > 0.0 for v in vals):
                raise InputError(f"{path}:{lineno}: mass but no valid radius")
            rows.append((float(r["Altitude"]), frac, *vals))
    alt = np.array([r[0] for r in rows])
    if np.any(np.diff(alt) <= 0.0):
        raise InputError(f"{path}: altitudes not strictly increasing")
    return [np.array(c) for c in zip(*rows)]


def read_rho_droplet(path):
    with open(path, newline="") as f:
        return {float(r["d0_um"]): float(r["rho_droplet_kg_m3"]) for r in csv.DictReader(f)}


def narrow_bins(r_nm):
    """r_edge and, per row, the index of its narrow bin (merging near-equal radii)."""
    order = np.argsort(r_nm)
    lnr = np.log(r_nm[order])
    group = np.concatenate([[0], np.cumsum(np.diff(lnr) > MERGE_LNR)])
    centres = np.exp(np.array([lnr[group == g].mean() for g in range(group[-1] + 1)]))
    edges = np.ravel(np.column_stack([centres * math.exp(-HALF_WIDTH),
                                      centres * math.exp(HALF_WIDTH)]))
    row_bin = np.empty(len(r_nm), dtype=int)
    row_bin[order] = 2 * group  # narrow bins are 0, 2, 4, ...; odd bins are empty gaps
    if np.any(np.diff(edges) <= 0.0) or np.max(np.abs(np.log(r_nm / centres[row_bin // 2]))) > MERGE_LNR:
        raise InputError("radii too close together to give each a narrow bin")
    return edges, row_bin, centres


def carma_placement(r_nm, weight, rmin_cm, vmrat, nbin):
    """Share of the weight per CARMA bin (1-based) and outside the grid."""
    _, rlow, rup = carma_bin_edges(rmin_cm, vmrat, nbin)
    b = np.searchsorted(rup, r_nm)
    outside = (r_nm <= rlow[0]) | (r_nm >= rup[-1])
    share = {int(i) + 1: weight[(b == i) & ~outside].sum() / weight.sum()
             for i in np.unique(b[~outside])}
    return share, weight[outside].sum() / weight.sum()


def write_nc(out, d0_um, alt, r_edge, row_bin, r_nm, r_um, rho_droplet, rho_alumina,
             radius_factor, max_spread, bin_width, source):
    nalt, nbin = len(alt), len(r_edge) - 1
    mass_frac = np.zeros((nalt, nbin))
    mass_frac[np.arange(nalt), row_bin] = 1.0
    with netcdf_file(out, "w", version=2) as nc:
        nc.createDimension("psd_alt", nalt)
        nc.createDimension("psd_bin", nbin)
        nc.createDimension("psd_bin_edge", nbin + 1)
        v = nc.createVariable("psd_alt", "f8", ("psd_alt",))
        v.units = "km"
        v.long_name = "altitude above ground of PSD rows"
        v[:] = alt
        v = nc.createVariable("r_edge", "f8", ("psd_bin_edge",))
        v.units = "nm"
        v.long_name = "PSD bin edge radius"
        v[:] = r_edge
        v = nc.createVariable("mass_frac", "f8", ("psd_alt", "psd_bin"))
        v.units = "1"
        v.long_name = "fraction of emitted mass in each PSD bin (rows sum to 1)"
        v[:] = mass_frac
        v = nc.createVariable("particle_radius", "f8", ("psd_alt",))
        v.units = "nm"
        v.long_name = "alumina particle radius of each row (mass-averaged at that altitude)"
        v[:] = r_nm
        v = nc.createVariable("droplet_radius", "f8", ("psd_alt",))
        v.units = "um"
        v.long_name = "mass-averaged final Al droplet radius at that altitude"
        v[:] = r_um
        nc.title = (f"CARMA alumina emission particle size distribution (monodisperse at each "
                    f"altitude), initial droplet diameter {d0_um:g} um")
        nc.input_size = f"{d0_um:g} um initial Al droplet diameter"
        nc.max_spread_ln_r = np.float64(max_spread)
        nc.carma_bin_width_ln_r = np.float64(bin_width)
        nc.rho_droplet_kg_m3 = np.float64(rho_droplet)
        nc.rho_alumina_kg_m3 = np.float64(rho_alumina)
        nc.al_to_al2o3_radius_factor = np.float64(radius_factor)
        nc.source_radius_profile = str(source)
        nc.history = "created by droplet_summary_to_alumina_psd.py"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile-dir", type=Path, default=PROFILE_DIR)
    p.add_argument("-o", "--outdir", type=Path, default=HERE / "droplet_size_psd")
    p.add_argument("--rho-alumina", type=float, default=3.95,
                   help="alumina particle density in g/cm3 (default 3.95, the density "
                        "in the CARMA alumina model)")
    p.add_argument("--vmrat", type=float, default=2.0,
                   help="CARMA bin volume ratio for the spread test (default 2)")
    args = p.parse_args(argv)

    rho_alumina = args.rho_alumina * 1000.0
    bin_width = math.log(args.vmrat) / 3.0
    rho_droplet = read_rho_droplet(args.profile_dir / "droplet_emission_summary.csv")
    sizes = sorted((float(m.group(1)), path) for path in args.profile_dir.iterdir()
                   if (m := RADIUS_RE.match(path.name)))
    if not sizes:
        raise InputError(f"{args.profile_dir}: no droplet_<D>um_final_radius_profile.csv")
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"spread test per altitude: ln(r_max / r_min) < ln(vmrat)/3 = {bin_width:.4f}")
    for d0, path in sizes:
        alt, frac, r_mean, r_min, r_max = read_radius_profile(path)
        if np.any(r_mean < r_min * (1 - RADIUS_RTOL)) or np.any(r_mean > r_max * (1 + RADIUS_RTOL)):
            raise InputError(f"{path}: mass-averaged radius outside [r_min, r_max]")
        spread = np.log(r_max / r_min)
        if spread.max() >= bin_width:
            a = alt[spread.argmax()]
            raise InputError(f"d0 = {d0:g} um, {a:g} km: radius spread {spread.max():.4f} in "
                             f"ln r is not narrower than one CARMA bin ({bin_width:.4f}); a "
                             "polydisperse row is needed")

        radius_factor = (rho_droplet[d0] / (AL_PER_AL2O3 * rho_alumina)) ** (1.0 / 3.0)
        r_nm = r_mean * radius_factor * UM_TO_NM
        r_edge, row_bin, centres = narrow_bins(r_nm)
        out = args.outdir / f"alumina_droplet_{size_tag(d0)}_psd.nc"
        write_nc(out, d0, alt, r_edge, row_bin, r_nm, r_mean, rho_droplet[d0], rho_alumina,
                 radius_factor, spread.max(), bin_width, path)

        r0_nm = d0 / 2.0 * UM_TO_NM
        mean_ratio = np.sum(frac * r_nm) / np.sum(frac) / r0_nm
        print(f"\nd0 = {d0:g} um -> {out.name}: {len(alt)} altitudes "
              f"({alt[0]:g}-{alt[-1]:g} km), {len(centres)} distinct radii, "
              f"max spread {spread.max():.4f}")
        print(f"  alumina radius {r_nm.min() / r0_nm:.4f}-{r_nm.max() / r0_nm:.4f} x initial "
              f"radius (emitted-mass mean {mean_ratio:.4f})")
        rmin_planned_nm = r0_nm * PLANNED_RMIN_FACTOR
        grids = [(f"alumina_{size_tag(d0)} grid (rmin {r0_nm:g} nm, {SIZE_MODEL_NBIN} bins)",
                  r0_nm * 1e-7, SIZE_MODEL_NBIN),
                 (f"planned alumina_{size_tag(d0)} grid (rmin {rmin_planned_nm:.6g} nm, "
                  f"{SIZE_MODEL_NBIN} bins)", rmin_planned_nm * 1e-7, SIZE_MODEL_NBIN),
                 (f"base alumina grid (rmin 0.5 nm, {BASE_NBIN} bins)", BASE_RMIN_CM, BASE_NBIN)]
        for name, rmin_cm, nbin in grids:
            share, clipped = carma_placement(r_nm, frac, rmin_cm, args.vmrat, nbin)
            if not share:
                print(f"  {name}: outside the grid")
                continue
            parts = ", ".join(f"bin {b} {s:.1%}" for b, s in share.items())
            print(f"  {name}: {parts}; outside the grid {clipped:.1%} of emitted mass")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except InputError as e:
        sys.exit(f"error: {e}")
