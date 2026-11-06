# Plan: lower rmin by one bin in the alumina_d<D>um size models

Written 2026-10-24 on the laptop, where the CAM tree is not available. Every CAM-side step below is
meant to be run on the machine that has the CAM tree, on branch `ashajain/carma_aluminum_model`,
with case `FWsc2010climo.CARMA_dev`. File names inside CAM that were not visible from the laptop
are found with `grep` rather than assumed; check each hit before editing.

## Why

The droplet PSD files (`droplet_size_psd/alumina_droplet_d<D>um_psd.nc`) give, at each emission
altitude, one alumina particle radius: the mass-averaged final radius of the droplets that stop
there, converted from Al (2.4 g/cm3) to Al2O3 (3.95 g/cm3). Those radii run from about 0.86 to 1.03
times the starting droplet radius r0.

The size models were built with `rmin = r0` (16 bins, vmrat 2). The lower edge of bin 1 is then
`r0 * (2/3)^(1/3) = 0.874 r0`, so 20 % to 53 % of the emitted mass (depending on the size) has a
radius below the grid. The PSD remap clips it; more than 1 % clipped in a row stops the model, and
a row that is entirely below the grid has nothing left to renormalise.

Lowering `rmin` by one bin, to `r0 * 2^(-1/3) = 0.794 r0`, moves the lower edge of bin 1 to
`0.693 r0`. Every emitted radius then lands in bin 1 or bin 2 and nothing is clipped
(`droplet_summary_to_alumina_psd.py` prints this as the "planned" grid).

Side effects:

- The old starting radius r0 becomes exactly the centre of the new bin 2. The existing
  `alumina_mono_d<D>um_psd.nc` files therefore still put all their mass in one bin at the same
  radius as before (bin 2 instead of bin 1), which gives a regression test (step 7).
- The number of bins stays 16, so PCNST (57), the CRALUM01-16 history fields and N_RAD_CNST
  should not change. The top of the grid moves down one bin, from about 35 r0 to about 28 r0
  (for example 3,494 um instead of 4,402 um for d250um). Emitted particles only shrink or stay the
  same, so the top bins should stay empty; if coagulation is on and you want the old top, use
  17 bins instead, but then the history field lists, the optics count and PCNST all change.

## New rmin values

All radii in cm. "Current rmin" is r0 = d/2, as recorded in CASE_LOG.md (2026-09-30).

| Model | Current rmin (cm) | New rmin (cm) | New bin 1 lower edge (nm) | New bin 2 centre (nm) | New bin 16 upper edge (um) |
|---|---|---|---|---|---|
| alumina_d2um | 1.0000e-04 | 7.937005260e-05 | 693.361 | 1000 | 27.955 |
| alumina_d6um | 3.0000e-04 | 2.381101578e-04 | 2080.08 | 3000 | 83.864 |
| alumina_d10um | 5.0000e-04 | 3.968502630e-04 | 3466.81 | 5000 | 139.77 |
| alumina_d14um | 7.0000e-04 | 5.555903682e-04 | 4853.53 | 7000 | 195.68 |
| alumina_d18um | 9.0000e-04 | 7.143304734e-04 | 6240.25 | 9000 | 251.59 |
| alumina_d25um (new model) | - | 9.921256575e-04 | 8667.02 | 12500 | 349.43 |
| alumina_d35um (new model) | - | 1.388975920e-03 | 12133.8 | 17500 | 489.21 |
| alumina_d45um (new model) | - | 1.785826183e-03 | 15600.6 | 22500 | 628.98 |
| alumina_d50um | 2.5000e-03 | 1.984251315e-03 | 17334 | 25000 | 698.86 |
| alumina_d100um | 5.0000e-03 | 3.968502630e-03 | 34668.1 | 50000 | 1397.7 |
| alumina_d150um | 7.5000e-03 | 5.952753945e-03 | 52002.1 | 75000 | 2096.6 |
| alumina_d200um | 1.0000e-02 | 7.937005260e-03 | 69336.1 | 100000 | 2795.5 |
| alumina_d250um | 1.2500e-02 | 9.921256575e-03 | 86670.2 | 125000 | 3494.3 |

The 0.01 and 0.1 um cases fit in the base 36-bin `alumina` model (rmin 5e-8 cm), with nothing
clipped and no change needed. The 2 um case fits there too, so the d2um model is optional.

## Steps on the CAM machine

Set these for the commands below (adjust to your paths):

```bash
CAM=~/MIT/cesm2_port/components/cam
CASE=~/MIT/cesm2_port/cases/FWsc2010climo.CARMA_dev
EMIS=~/MIT/cesm2_port/emissions
```

### 0. Check the PSD reader can take the new file sizes

The droplet PSD files have 95 to 121 altitude rows and 37 to 83 PSD bins. The only PSD file run so
far (`alumina_test_psd.nc`) had 9 rows and 20 bins. Look for fixed-size arrays or limits:

```bash
grep -n -i -E "parameter|allocate|dimension|max_|psd_alt|psd_bin" $(grep -rl "alumina_psd_binfrac" $CAM/src)
```

If the arrays are allocated from the file's dimensions, nothing to do. If there is a fixed limit
below 121 rows or 83 bins, raise it first.

### 1. Find every place the size models are defined

```bash
grep -rIl "alumina_d250um" $CAM | sort -u
grep -n -i "rmin" $CAM/src/physics/carma/models/alumina_d250um/carma_model_mod.F90
```

Expected (from CASE_LOG.md): one `carma_model_mod.F90` per model holding rmin; the other modules
and `vertical.F90` are symlinks to `../alumina`; build-namelist sets `carma_do_optics = .false.`,
the pre-built optics paths and the default PSD file for each model. Note every file the first
grep lists, since the new 25, 35 and 45 um models need an entry in each.

### 2. Change rmin in the 10 existing models

In each `models/alumina_d<D>um/carma_model_mod.F90`, replace the rmin value with the "New rmin"
from the table. Write it as a literal, as the current files do, or as `r0 * 2._f**(-1._f/3._f)`
if the file already computes rmin from the radius. After editing, confirm:

```bash
for d in 2 6 10 14 18 50 100 150 200 250; do
  echo -n "d${d}um: "; grep -n -i "rmin *=" $CAM/src/physics/carma/models/alumina_d${d}um/carma_model_mod.F90
done
```

### 3. Create alumina_d25um, alumina_d35um and alumina_d45um

Copy the d50um directory the same way the others were made (only `carma_model_mod.F90` is a real
copy; keep the symlinks relative):

```bash
cd $CAM/src/physics/carma/models
for d in 25 35 45; do
  mkdir alumina_d${d}um
  cp alumina_d50um/carma_model_mod.F90 alumina_d${d}um/
  for f in $(find alumina_d50um -type l -exec basename {} \;); do
    ln -s "$(readlink alumina_d50um/$f)" alumina_d${d}um/$f
  done
done
ls -l alumina_d25um alumina_d35um alumina_d45um
```

Set each new model's rmin from the table, and change any other `d50um` text in the copied
`carma_model_mod.F90` (model name, comments). Then add the three models to every file found in
step 1 (build-namelist defaults, optics paths, default PSD file), copying the d50um entries.

### 4. Regenerate the pre-built optics

The bin radii change, so all optics files for the size models must be rebuilt: 13 models x 16 bins
= 208 files (the old 160 plus 48 for the new models). Move the old files aside instead of
deleting them:

```bash
OPT=<inputdata>/atm/cam/physprops/alumina
mkdir -p $OPT/rmin_r0_2026-09-30
mv $OPT/alumina_d*um_CRALUM* $OPT/rmin_r0_2026-09-30/
```

Then run `$CAM/tools/carma_optics_offline` for each of the 13 models the same way as on
2026-09-30 (see the tool's README for how it takes the model or rmin). Check that it picked up the
new rmin: the bin-1 radius in the d250um files should be 99,212.6 nm (was 125,000 nm), and the
count should be 208. As before, the patched bhmie is needed for the large-size, |m| < 1 bands.

### 5. Put the new emission and PSD files where the model can read them

From this repo (`alumina_emissions`), copy `droplet_size_events/*.nc` and `droplet_size_psd/*.nc`
to `<inputdata>/atm/waccm/emis/` (or `$EMIS`), and check the PSD file with the same remap as the
Fortran:

```bash
python remap_psd_reference.py droplet_size_psd/alumina_droplet_d250um_psd.nc --rmin 9.921256575e-03 --nbin 16
```

It should report a clipped fraction of 0 at every altitude, with mass only in bins 1 and 2.
(`remap_psd_reference.py` needs netCDF4, which the WACCM env has.)

### 6. Build and preview

For one model first (d250um):

```bash
cd $CASE
./xmlchange CAM_CONFIG_OPTS="-phys cam6 -age_of_air_trcs -chem waccm_sc_mam4 -carma alumina_d250um"
./preview_namelists
./case.build --clean atm && ./case.build
```

Expect PCNST = 57 and N_RAD_CNST = 64, as before; check them wherever you read them on
2026-09-30 (the preview_namelists/build output). Repeat the preview for the other 12 models
(no run needed) to check each finds its optics and PSD files.

### 7. Regression run: old files on the new grid

With the old emission file and the old monodisperse PSD, the particles have exactly the same radius
as before (r0, now the centre of bin 2), so the 2-day d250um run of 2026-10-01 should be
reproduced:

- `carma_emis_file` = `alumina_2025_reentry_events.nc`
- `carma_emis_psd_file` = `alumina_mono_d250um_psd.nc`
- 2 days, `carma_emis_scale` = 1.0

Expected: emitted 1,897.950869 kg (day 1: 7 events, 1,781.95 kg; day 2: 1 event, 116.00 kg),
all of it in bin 2 (bin 1 empty); DD 1,852.15 kg and SW 45.80 kg; end-of-day-2 burden 0. Any
difference beyond round-off means something other than rmin changed.

### 8. Test run: droplet files on the new grid

- `carma_emis_file` = `alumina_2025_reentry_events_d250um.nc`
- `carma_emis_psd_file` = `alumina_droplet_d250um_psd.nc`
- 2 days, `carma_emis_psd_allow_clip` left false

Expected emitted mass (from the files):

| Case | Day 1 (7 events) | Day 2 (1 event) | 2-day total |
|---|---|---|---|
| d250um | 1,211.274672 kg: bin 1 239.6398, bin 2 971.6349 | 78.850620 kg: bin 1 15.5999, bin 2 63.2507 | 1,290.125292 kg |
| d10um | 1,187.386471 kg: bin 1 419.0764, bin 2 768.3100 | 77.295564 kg: bin 1 27.2807, bin 2 50.0149 | 1,264.682035 kg |

The d250um total is the old file's 1,897.950872 kg times the d250um surviving share (0.679747).
No endrun from clipping, and no mass in bins 3-16 at emission.

### 9. Commit on the CAM branch

Commit the rmin changes, the three new model directories and the build-namelist entries on
`ashajain/carma_aluminum_model`, and record the runs in CASE_LOG.md here.
