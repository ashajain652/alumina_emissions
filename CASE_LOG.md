# Change log for case FWsc2010climo.CARMA_dev

Original state (2026-09-27, before any CARMA work):
- CAM_CONFIG_OPTS = `-phys cam6 -age_of_air_trcs -chem waccm_sc_mam4`
- user_nl_cam: `nhtfrq = -24`, `mfilt = 1`
- STOP_OPTION=ndays, STOP_N=2, NTASKS_ATM=12, RUN_STARTDATE=0001-01-01
- The last successful build and run (st_archive success) was at 2026-09-27 11:28.

## WP0: meteor_smoke baseline
- `./xmlchange --append CAM_CONFIG_OPTS=" -carma meteor_smoke"`
- `./case.build --clean atm && ./case.build`
- Downloaded meteor_smoke_kalashnikova.nc and smoke_grf_frentzke.nc with `./check_input_data --download`
- Build succeeded (models/meteor_smoke in Filepath)
- `./case.submit` (2 days)
- Run completed at 15:24 in about 8 min wall time, and st_archive succeeded. DUST01–DUST28 are present in h0.
- Check: the daily-mean DUST burdens are 22,374 kg (day 1) and 66,210 kg (day 2), against 21,918 and 65,753 kg expected for 16 kt/yr (burden growing linearly from 0). They agree to within 2 %. **WP0 passed.**

## WP5: alumina integration
- `./xmlchange CAM_CONFIG_OPTS="-phys cam6 -age_of_air_trcs -chem waccm_sc_mam4 -carma alumina"`
- user_nl_cam rewritten: the original nhtfrq=-24 and mfilt=1 are kept for h0; added the carma_emis_file and carma_emis_psd_file paths, plus fincl2 (daily-mean EM/DD/SW) and fincl3 (end-of-day instantaneous CRALUM, PS).
- Build succeeded (PCNST=77, no alumina warnings). `./case.submit` (2 days)
- Run completed at 15:36 in about 8 min, and st_archive succeeded.
- **WP6 results (2-day test, 10 events, 90 rows, 1693.9828 kg):**
  1. Firing: all 90 rows fired exactly once, each in the step whose window [prev, curr) contains its time. No endrun.
  2. Mass: emitted (daily-mean EM with h0 PS) = 1693.982797 kg, relative error 2e-9. End-of-day burden from h2: day 1 = 882.2233 kg (exactly events 1–5); day 2 = 1693.982800 kg. DD+SW deposition = 0.
  3. PSD remap: the Fortran psdfrac matches remap_psd_reference.py to all printed digits. Clipped fraction 9e-16.
  4. EM size split at all 73 firing cells matches the reference to within 1.4e-8.
  - Mapping: every row was placed within 113 km of its location. No NaN or negative values. 99.6 % of the mass is above 1 hPa and none is below 10 hPa (mmr there ≤ 1e-24, numerical noise).
  - Not run yet: check 5 (zero source) and check 6 (a full rebuild and run of meteor_smoke with the new bld/ code).
- `'PS'` added to fincl2, so emitted_mass.py works on h1 directly in future runs.

## v2: alumina radiatively active (2026-09-30)
- Refractive indices: `emissions/make_alumina_refidx.py` interpolates n (tropf_real_alumina_2.csv) and k (tropf_kapoor_combined_img_alumina.csv) to the 30 CARMA band centres (band-centre values, no band averaging). The table is hardcoded in models/alumina/carma_model_mod.F90; see alumina_refidx_bands.csv.
- New namelist switch `carma_alumina_rad_feedback` (default .true.): CRALUM01-36 go into rad_climate, and rad_diag_1 = the same list without alumina. If .false., the reverse.
- configure raises max_n_rad_cnst from 30 to 64 for alumina (rad_climate now has 51 entries).
- user_nl_cam: added a 4th tape (h3, daily mean): FSNT/FLNT/FSNS/FLNS/FSUTOA/FLUT/QRS/QRL and their _d1 versions, AEROD_v, AODvstrt, ODV_CRALUM01-36. The previous user_nl_cam was saved in the session scratchpad.
- Build (clean atm): OK, N_RAD_CNST=64. Run 1 (2 days): OK in about 10 min. CARMA wrote alumina_CRALUM01-36_rrtmg.nc; Postprocessing/alumina/check_alumina_optics.py matches an offline Mie calculation to within 4e-7 (refidx exact).
- **Problem found in run 1:** X - X_d1 (climate with alumina minus diagnostic without) gave -0.012 W/m2 global mean and up to -1.35 W/m2 locally, with the largest heating differences near 990 hPa. That is roughly 5 orders of magnitude more than alumina with AOD <= 1.3e-6 can produce. Cause: modal_aer_opt.F90 recomputes MAM size and water uptake for diagnostic lists (list_idx > 0), while the climate call uses the pbuf values, so the climate and diagnostic calls differ even with identical MAM entries.
- **Fix:** build-namelist now sets rad_diag_1 = the list with alumina and rad_diag_2 = the list without, and rad_climate includes alumina only if carma_alumina_rad_feedback is true. The alumina effect is X_d1 - X_d2. fincl4 now also has the _d2 fields.
- Run 2 stopped at startup to add ndens = 2,2,2,1 (h3 in double precision): the expected alumina effect (~1e-5 W/m2) is below float32 resolution of fluxes (~3e-5 at 300 W/m2).
- Run 2 (2 days, rad_diag_1/2 pair, double-precision h3): OK at 14:44, and st_archive succeeded. Results (`Postprocessing/alumina/alumina_rad_effect.py`):
  - Null check: in the first h3 record (no alumina yet), d1 - d2 = 0 exactly in every field. The old X - X_d1 artifact is still present (-0.012 W/m2), which confirms the diagnosis.
  - Day 1 / day 2 global mean dFSNT = -1.27e-7 / -3.92e-7 W/m2; dFLNT = -5.3e-9 / -9.8e-9 W/m2.
  - Largest local dFSNT = -1.1e-4 W/m2 at 71.1S 250E, the same cell as the largest alumina visible AOD (1.27e-6). The heating perturbations peak at 0.02-0.05 hPa (about 70 km), where the particles are.
  - Regression: emitted mass 1693.982797 kg (rel. err. -1.9e-9); end-of-day-2 burden 1693.982800 kg. Unchanged from v1.

## Particle size sensitivity models (2026-09-30)
- New CARMA models alumina_d{2,6,10,14,18,50,100,150,200,250}um: 16 bins, vmrat 2, rmin = emission radius (d/2). Only carma_model_mod.F90 is copied; the other three modules are relative symlinks to ../alumina. The default PSD file is atm/waccm/emis/alumina_mono_d<D>um_psd.nc; the 10 PSD files were copied from emissions/psd_sensitivitytest to inputdata/atm/waccm/emis/.
- Checks run in this case (temporary): user_nl_cam was set aside (a copy is in the session scratchpad) and CAM_CONFIG_OPTS switched to each -carma alumina_d<D>um for preview_namelists. All 10 gave PCNST=57, N_RAD_CNST=64, their own models/ dir and PSD file, and 16 CRALUM in rad_climate and rad_diag_1 (0 in rad_diag_2).
- Compile test: `./case.build --clean atm && ./case.build` with -carma alumina_d250um.
- Build OK: all four alumina_d250um sources, including the three symlinks, compiled from models/alumina_d250um/. The case was then restored: CAM_CONFIG_OPTS back to -carma alumina, original user_nl_cam, clean atm rebuild OK (PCNST=77, N_RAD_CNST=64).

## alumina_d250um short test (2026-09-30)
- The 36-bin user_nl_cam for the base alumina model is saved as emissions/user_nl_cam.alumina_base.
- user_nl_cam trimmed for 16-bin models: removed the carma_emis_psd_file override (the per-model default is used); fincl2/3/4 cut to CRALUM01-16.
- `./xmlchange CAM_CONFIG_OPTS="-phys cam6 -age_of_air_trcs -chem waccm_sc_mam4 -carma alumina_d250um",STOP_N=1`; clean atm rebuild; 1-day run.
- **Run 1 (d250um) failed at startup:** `miess:: The upper limit for acap is not enough`, in bin 3 at LW07. CARMA's default Mie routine caps its series at |m|x terms, and alumina has |m| < 1 in LW06 and LW07. CARMA wrote only bins 1-2 (bin 3 was an empty file).
- **Fix: pre-built optics.** tools/carma_optics_offline compiles CARMA's base/ code (with CAM's carma_precision and constants modules and small CAM stand-ins) into an offline driver.
  - Validated: offline miess is bit-identical to the 36 base-alumina files CARMA wrote at runtime.
  - CARMA's bhmie had an unconverged d(n) recurrence start (up to 0.3% error in Qsca). A patched copy in the tool (start at 2*max(xstop,|mx|)+15) matches miess to <= 5e-5 everywhere miess works: all bins of the base model and of d2/d6/d10, plus d250 bins 1-2.
  - 160 files (10 models x 16 bins) written to inputdata/atm/cam/physprops/alumina/.
  - build-namelist: the alumina_d<D>um models use carma_do_optics = .false. and absolute paths to those files. The partial alumina_d250um_CRALUM01-03 files were deleted from the run dir.
- Run 2 (d250um, 1 day): submitted.
- **Run 2 (d250um, 1 day): OK, but no sedimentation.** The pre-built optics were read from inputdata and CARMA wrote none. Mass was conserved (882.2233 kg emitted = airborne), but DD+SW = 0 and 99.8% of the mass was still above 1 hPa (55-80 km). Cause: base/vertical.F90 skips sedimentation for a column whose peak number concentration is at or below FEW_PC (CAM sets 1e-6 /cm3); the d250um peak was 5.9e-12 /cm3.
- **Fix:** models/alumina/vertical.F90 overrides CARMA base/vertical.F90 (the model dir comes first in Filepath) and drops the FEW_PC test for sedimentation only. The 10 size models symlink to it. base/ is an SVN external (carma3_49_rel), not in CAM git, so it is left unmodified (svn revert).
- Clean atm rebuild OK; the override was compiled from models/alumina_d250um/vertical.F90.
- user_nl_cam: carma_emis_scale = 1000 (user choice); 'CRALUMVM' (group fall velocity) added to fincl2.
- Run 3 (d250um, 1 day, x1000): submitted.
- Commits (CAM branch ashajain/carma_aluminum_model, not pushed): d228bdf2 (build-namelist uses pre-built optics for the size models) and 60d657f4 (offline optics tool). The tool moved from ~/MIT/cesm2_port/tools/ to components/cam/tools/carma_optics_offline (paths now relative to the CAM tree). Regression after the move: identical output. The models/alumina/vertical.F90 override is not committed yet, pending run 3.
- Run 3 result: sedimentation works (86% of 882 t deposited in 1 day: DD 742.9 t, SW 15.2 t, airborne 159.0 t). But the budget is off by +34.8 t (+3.9%), all in bin 1. Run 4: the same, plus a temporary h4 of every-step bin-1 budget terms to find the source (the run-3 user_nl_cam is saved as emissions/user_nl_cam.d250um_run3).
- **Run 4 (d250um x1000, every-step bin-1 budget in h4): the 3.9% is a timing artifact, not a loss of mass.** Instantaneous snapshots are taken before CARMA's tendencies, so the end-of-day burden misses exactly the last step's DD (34.8373 t). Every-step budget closes to 5e-15 relative. CRALUMVM (24-2000 m/s) is contaminated by the SMALL_PC floor in empty large bins; CARMA's bin-1 speed is 2.5 m/s near the surface and 121 m/s at 70 km (setupvf_std). Notes added to Postprocessing/alumina/README.md.
- Commit b50a5d74 (vertical.F90 override). Temporary h4 removed from user_nl_cam (back to the run-3 version: d250um settings, carma_emis_scale = 1000, CRALUMVM in fincl2). STOP_N set back to 2. Case left on -carma alumina_d250um (d250um build). The h4 files (0001-01-01 and -02, about 775 MB) remain in the archive.

## alumina_d250um with 2025 reentry emissions (2026-10-01)
- New emission file: emissions/alumina_2025_reentry_events.nc. It is built by reentry_csv_to_alumina_emissions.py, which applies ESA ablation profiles integrated over 5 km bands from 30 to 90 km, and csv_to_alumina_events.py --year-offset -2024. It has 1223 located reentries and 192,194.64 kg in total.
- user_nl_cam: carma_emis_file now points to that file, and carma_emis_scale = 1000 is removed (so the default 1.0 applies). Everything else is unchanged from user_nl_cam.d250um_run3. Run 4's archived history files were moved to archive/.../atm/hist/run4_d250um_x1000/.
- Run (d250um, 2 days, no rebuild): OK at 12:07, and st_archive succeeded. The model read 15899 rows with carma_emis_scale = 1.0.
  - Emitted: 1897.950869 kg (day 1: 7 events, 1781.95 kg; day 2: 1 event, 116.00 kg). Expected 1897.950872 kg, relative error -1.4e-9.
  - Deposited: DD 1852.15 kg + SW 45.80 kg = 1897.950906 kg. The end-of-day-2 burden is 0, so the 250 um particles all reach the surface within the 2 days. The end-of-day-1 snapshot holds 9.23 kg.
- Case left on -carma alumina_d250um (d250um build) with the 2025 user_nl_cam, STOP_N=2. Run output: archive/.../atm/hist/*.cam.h[0-3,6].0001-01-0{1,2,3}; submit log: emissions/d250um_2025_run.log.

## Fall velocity table test: FWsc2010climo.alum2um_vftest (2026-10-24)
- New case, cloned from FWsc2010climo.CARMA_dev (new build): -carma alumina_d2um, CAM branch ashajain/aerosol_terminal_velocity (uncommitted carma_vf_from_file / carma_vf_file code), same 2025 user_nl_cam plus `carma_vf_from_file = .true.` and `carma_vf_file = emissions/vf_tables/alumina_d2um_vf_carma.nc`. 2 days.
- The table holds CARMA's own layer-centre vf for d2um on this case's reference atmosphere (offline driver from the CARMA source, T from carma_reft.nc), at the layer altitudes alumina_vf_apply computes, so the used vf should equal CARMA's.
- Build OK (setupvf.o from models/alumina_d2um, calls alumina_vf_apply). Run OK at 15:05, and st_archive succeeded.
  - atm.log: table read (70 altitudes, 0.055-142.889 km, 16 bins); alumina_vf_apply listed all 1120 (lev, bin) pairs. Used vs CARMA and used vs table both agree to 9e-5, the log's print precision; the logged layer altitudes equal the table's.
  - Emissions fired (104 row firings in cesm.log). End-of-day burden 1781.951 kg (day 1) and 1897.951 kg (day 2), equal to the emitted mass; 2 um particles barely sediment in 2 days. No NaN.
- Found: CAM's cam/carma_getT.F90 reads the IC T into init_t(plon,pver,plat) but the file is (lon,lat,lev), so carma_reft.nc (the fixedinit reference T) is a scrambled 252-260 K at every level instead of the real profile (e.g. 794 K at the top, ~193 K near the mesopause). Reproduced to 1e-12 by re-reading the IC in that order. Not fixed yet.

## Droplet-size emission files for 2025 reentries (2026-10-24)
- New script: droplet_profiles_to_alumina_events.py. It writes 15 event files, droplet_size_events/alumina_2025_reentry_events_d<D>um.nc, one per droplet diameter (0.01, 0.1, 2, 6, 10, 14, 18, 25, 35, 45, 50, 100, 150, 200, 250 um). Each file's title gives the assumed particle diameter. The variables and the format (64-bit offset) are the same as alumina_2025_reentry_events.nc. The files are written with scipy, because no local env has netCDF4.
- The profiles come from Reentry Emissions/Droplet Emission Profiles/run_droplet_emission_profiles.py (outside this repo). Updated today: the ESA satellite profile is normalised to 1 (its raw integral is 0.5918), so it gives only the release shape. Each 1 m release slice is multiplied by the droplet's surviving mass fraction (1 - mass lost to vaporisation) and moved to the droplet's final altitude, then binned per km. A profile therefore sums to the surviving share, 0.666 (0.01 um) to 0.680 (250 um).
- Emission per row = inventory Al2O3 (2025_emissions.csv, 1223 located reentries, 192,194.64 kg) x profile value. Altitudes are the 1 km bin centres, 0.5-120.5 km; zero bins are dropped. Rocket bodies (S0-S4, 78.4 t) use the satellite profile. Options: --satellites-only drops rocket bodies; --conserve-mass normalises each profile to 1 (no vaporisation loss).
- Totals run from 128,022.81 kg (0.01 um) to 130,643.68 kg (250 um). Each file has 116,185-147,983 rows (5-6.5 MB), against 15,899 rows in the 5 km band file. Not yet tested in CARMA.
- Checks: event ids, dates, times and lat/lon match alumina_2025_reentry_events.nc. Per-event mass = inventory x profile sum to 7e-15 relative. The files from the earlier version (division by the ESA integral done in this script) agree to 7e-10.
- Only 10 of the 15 sizes have alumina_d<D>um CARMA models; 0.01, 0.1, 25, 35 and 45 um have PSD files only.
- Al2O3/Al conversion: the script now converts inventory Al2O3 to Al (x 2 M_Al / M_Al2O3 = 0.529257, AA7075 droplets treated as pure Al), applies the droplet profile to the Al, and converts each bin back to Al2O3. 192,194.64 kg Al2O3 = 101,720.35 kg Al. The profile is a mass fraction, so the conversions cancel: the regenerated files match the previous ones to 3e-16 relative (totals unchanged). The new global attribute al_al2o3_conversion records this.

## Droplet-size PSD files and size-model rmin plan (2026-10-24)
- New script: droplet_summary_to_alumina_psd.py writes droplet_size_psd/alumina_droplet_d<D>um_psd.nc (15 sizes). Each emission altitude gets one row, monodisperse at the mass-averaged final droplet radius there, converted from Al (2.4 g/cm3) to Al2O3 (3.95 g/cm3, the CARMA alumina density). psd_alt equals the emission files' altitudes exactly, so the model's altitude interpolation never blends rows. Input: Droplet Emission Profiles/droplet_<D>um_final_radius_profile.csv from run_droplet_emission_profiles.py (outside this repo).
- Per-altitude radius spread is at most 0.087 in ln r (one CARMA bin is 0.231). Alumina radii run from 0.86 to 1.03 x the initial radius r0.
- Problem: the size models have rmin = r0, so bin 1 starts at 0.874 r0 and 20-53 % of the emitted mass would be clipped (endrun). The base 36-bin grid (0.01, 0.1, 2 um) has no clipping.
- Plan (SIZE_MODEL_RMIN_PLAN.md, not yet run): rmin = r0 x 2^(-1/3) in all size models, new d25/d35/d45 models, regenerate 208 optics files, a regression run (old files must reproduce 1897.950869 kg emitted, DD 1852.15 kg, SW 45.80 kg) and a droplet-file test (d250um: 1290.125292 kg over 2 days). With the planned grid nothing is clipped; mass splits between bins 1 and 2.
