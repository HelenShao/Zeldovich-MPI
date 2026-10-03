# Production neutrino ICs: FLAMINGO-phase-matched Abacus runs (cosm202-206)

Last updated: 2026-09-30

Goal: run abacus_cosm202-206 (the massive-neutrino cosmologies) with the same initial white noise as
FLAMINGO L1000N1800, each coloured with its own P(k), so the Abacus runs can be compared with FLAMINGO
object by object. The emulator-box runs (1215 Mpc/h, 4480^3) are deferred (section 6).

## 1. Why FLAMINGO matching requires FLAMINGO's box

Panphasia white noise is defined on integer mode indices n = (n_x, n_y, n_z) and has no box
dependence: `PANPHASIA_compute_kspace_field_` takes only the level and grid size. A run turns mode n
into the physical wavenumber k = 2 pi n / L. So the same noise in a different box puts the same
phases on different physical scales, and the structures no longer correspond. Phase matching to
FLAMINGO therefore means FLAMINGO's box (1000 Mpc) and its descriptor, with any N up to N0. The
emulator boxes (1215 Mpc/h) cannot be matched to FLAMINGO this way.

Because the noise is box independent, one set of noise files serves all five cosmologies, including
cosm206 with its different h and box.

## 2. Box and h convention

FLAMINGO's box is L = 1000 Mpc (not Mpc/h). Abacus works in Mpc/h, so `BoxSize = 1000 h`:

| Cosmology | h | BoxSize [Mpc/h] | Omega_Smooth | ZD_f_cluster = 1 - Omega_Smooth/Omega_M |
|---|---|---|---|---|
| cosm202 | 0.673 | 673.0 | 0.001422 | 0.99550 |
| cosm203 | 0.673 | 673.0 | 0.002845 | 0.99100 |
| cosm204 | 0.673 | 673.0 | 0.005689 | 0.98200 |
| cosm205 | 0.673 | 673.0 | 0.011378 | 0.96399 |
| cosm206 | 0.662 | 662.0 | 0.005880 | 0.98207 |

h is each cosmology's grid value (decided 2026-09-30). FLAMINGO's own h = 0.6732
(`FLAMINGO_test.conf`, `BoxLength = 673.2`) is not used, so our boxes differ from FLAMINGO's by
0.03 % in Mpc/h (cosm202-205). cosm206 corresponds to FLAMINGO's PlanckNu0p24Var box (h = 0.662).

## 3. Pipeline

1. **Noise, once.** monofonIC (Willem Elbers' `output_k_space_slabs` branch) with the FLAMINGO
   descriptor `[Panph6,L18,(56034,71400,250000),S1,KK1025,CH-999,Flamino_Gpc1]` and `GridRes = 1800`.
   Panphasia picks relative level 10, so N0 = 2 * 2^10 = 2048. Output: 1025 slabs
   `output_k_space_slab_y.{0..1024}`, 64 MiB each (about 69 GB), in
   `/flare/Abacus/helenshao/ICs/flamingo_panphasia_N2048/noise/` (read-only after the run, with
   `noise_sha256.txt`). Config: `flamingo_noise/FLAMINGO_noise_N1800.conf` (`LPTorder = 1`, test
   `white_noise`, so no particle ICs); job `flamingo_noise/pbs_flamingo_noise.sh`.
   - Caveat: the slab writer's index is a 32-bit `int` equal to `iy * N0 * (N0/2+1) + ...`. It overflows
     if one rank holds more than about 1000 y-planes at N0 = 2048, so always use at least 3 MPI ranks.
   - Run 2026-09-30 (job 8881605, 2 nodes, 24 ranks). All 1025 slabs were written at the correct size.
     monofonIC then aborted in its testing-mode `output_white_noise`, because of an off-by-one in
     `Grid_FFT::Compute_PowerSpectrum`: `nbins = int(kmax/kmin)` rounds to 899 at GridRes 1800, but the
     bin index reaches 899, so it writes one bin past the end and corrupts the heap. This happens after
     the slab write, and the writer only reads `return_field`.
     - The slabs are intact. Their low-k modes equal those of the same descriptor at N0 = 512 (the
       crop300 noise): r(k) = 1.00000 and amplitude ratio 1.00000 for every k <= 64. This is Panphasia's
       multi-resolution property, which also holds between N0 = 256 and 512.
     - `mono_N1800_white_noise.hdf5` was not written, so crop (b) skips its real-space block.
       Crop (a) already compares with monofonIC's own crop.
     - Post-steps (sha256, read-only) were run by hand; the core files are in `crash_8881605_cores/`.
2. **P(k).** CLASS v2.9.0 P_cb at z = 31 (and z = 49 as a check) for each cosmology, in
   `AbacusAurora/Cosmologies/FLAMINGO_cosm/abacus_cosm20N/`. The ini is the grid's
   `abacus_cosm20N.parameters.ini` with only `root` and `z_pk` changed. Job `run_class_z31.pbs`;
   `check_z49.py` writes `z49_check.txt`. The run of 2026-09-30 passes: all five z = 49 P_cb are
   bitwise identical to the grid's z = 49 output, and cosm202 z = 31 is identical to the earlier
   validated `pk/cosm202_z31`. Each subdirectory's `cosm.def` is a copy of the grid one
   with `ZD_Pk_filename` pointing at `abacus_cosm20N_z31.z1_pk_cb.dat` and `ZD_Pk_file_redshift = 31.0`.
3. **ICs, embedded.** A custom `multistep`: abacus `9575f04` (the production commit of the AbacusAurora
   runs) with `subprojects/zeldovich_mpi` replaced by branch `neutrinos`, built with
   `-Dzeldovich_mpi:extra_cpp_args=['-DLOAD_EXTERNAL_NOISE=1','-DEXTERNAL_NOISE_DIR="<noise dir>"']`
   on the 26.181.0 image, otherwise the production options. Since 2026-09-30, 26.181.0 is the default
   compute image; `next-eval` has a single node, so the jobs use `capacity` / `debug` / `debug-scaling`.
   Location:
   `/home/helenshao/abacus-store/9575f04_zmpi-neutrinos` (`custom_build/prepare_checkout.sh` on
   login, then `custom_build/pbs_build_multistep_nu.sh`). `provenance/` records the zeldovich_mpi
   commit and the uncommitted diff.
   - The noise directory is compiled in, and the environment variable `ZD_EXTERNAL_NOISE_DIR` (if set)
     overrides it, so the same binary serves the rockV2 round trip.
   - Zeldovich-MPI reads each y-plane's slab and keeps the modes with |k_i| < 900 by signed k (index
     k + 2048 for k < 0). This is exactly the mode set of a native 1800 grid; Nyquist planes are zeroed.
     D(k) = sqrt(P(k)) w(k).
4. **Runs.** `AbacusAurora/Simulations/FLAMINGO_nu_c20N.par2`: `NP = 1800**3`, `CPD = 441`
   (PPC 68.6), `NumZRanks = 2` on 2 nodes, z_init = 31, no light cones, `FinalRedshift = 0.0`.
   All five parse, and Abacus's `setup_zeldovich_params` gives `ZD_Pk_sigma_ratio = 1.0` exactly and
   accepts the `ZD_f_cluster` values (`flamingo_par2/setup_zeldovich_check.txt`).
   - The production binary would silently make seed-based ICs from `ZD_Seed = 0`. Check that the log
     shows `[external noise] reading .../flamingo_panphasia_N2048/noise/...`.

## 4. Physics choices

| Choice | Decision | Notes |
|---|---|---|
| Start redshift | z_init = 31 | FLAMINGO's start; matches the validated `pk/cosm202_z31.*`. |
| P(k) normalisation | `ZD_Pk_norm = 0` | The z = 31 file is used as-is. With `ZD_Pk_file_redshift = InitialRedshift = 31`, Abacus sets `ZD_Pk_sigma_ratio` = D(31)/D(31) = 1, so even `Pk_norm = 8` would give factor 1 (`power_spectrum.cpp`, lines 210-233); setting 0 removes the dependence on that coincidence. |
| P_cb(z_init) | CLASS P_cb at z = 31 | FLAMINGO's monofonIC used zwindstroom: z = 0 spectra back-scaled with scale-dependent 3-fluid growth, in N-body gauge. Open: the two differ at large scales (gauge, radiation, neutrino clustering; about 0.4 % variation of the growth between k_min and k = 0.1 at z = 31). This matters if the z = 0 amplitudes must match FLAMINGO to better than a percent. |
| LPT order | Abacus 2LPT (`LagrangianPTOrder = 2`) | FLAMINGO used 3LPT. Abacus's 3LPT asserts `Omega_Smooth == 0` (`lpt.cpp`), so it is unavailable for neutrino cosmologies. The 3LPT term scales as D^3 and is small at z = 31; the difference shows up as small early transients. |
| PLT | `ZD_qPLT = 1`, `ZD_qPLT_rescale = 1`, `ZD_PLT_target_z = 12` | The Abacus standard, as in the emulators (decided 2026-09-30). It changes only modes near Nyquist: eigenvector direction, per-mode growth rate f = (sqrt(1 + 24 e.val f_cluster) - 1)/4, and the rescale (a_12/a_31)^(f_fluid - f_PLT). Low-k modes stay exactly sqrt(P) w, so the phase match holds. Near Nyquist (k_Nyq ≈ 8.4 h/Mpc) our ICs differ from FLAMINGO's by design. |
| Neutrinos | Smooth (`Omega_Smooth`), `f_cluster = 1 - Omega_Smooth/Omega_M` | Abacus does not simulate neutrino particles. |
| Corner modes | `ZD_CornerModes = 1` with `ZD_k_cutoff = 1` | Decided 2026-10-01. The default (0) zeroes modes with |k| >= N/2 (the cube corners), which monofonIC keeps; 1 keeps the full cube of |k_i| < N/2 so the mode set matches FLAMINGO's. The emulators use the default. |
| Outputs | No light cones; `FinalRedshift = 0.0` | z = 0 for FLAMINGO snapshot comparisons (`Small.par2` stops at 0.1). |

## 5. Comparing with FLAMINGO

- Axes: Abacus (x, y, z) = monofonIC (z, y, x). Swap x and z when matching positions or fields.
- Sign: monofonIC's testing-mode `delta` has the opposite sign to our density (`--mono-sign -1` in
  `validate_external_noise.py`). Check the sign convention of FLAMINGO's particle ICs against ours
  (e.g. correlate early density fields) before object-by-object comparisons.
- Box: 673.0 vs 673.2 Mpc/h (h = 0.673 vs 0.6732) for cosm202-205.
- Redshifts: FLAMINGO keeps full particle snapshots at z = 5, 4, 3, 2, 1.5, 1, 0.75, 0.5, 0.4, 0.3, 0.2,
  0.1, 0 (all 78 for L1_m9 and L1_m9_DMO), halo catalogues at all 78 outputs (dz = 0.05 for z <= 3, 0.25
  for 3-5) and power spectra at 123 redshifts. Our subsamples (base.def: 0.1, 0.3, 0.5, 0.8, 1.1, 1.4,
  2.0, 2.35, 3.0, 3.5, 4.0) match six full snapshots and all fall on the catalogue grid. No z = 0
  subsample (decided 2026-10-01); the z = 0 state is in each run's final checkpoint. Gravity-only
  counterparts (`*_DMO`) exist for every cosmology variation.

### Production runs (`production/pbs_flamingo_nu.sh`, `qsub -v COSM=20N`, 2 nodes, `capacity`)

| Cosmology | Job | Result |
|-----------|-----|--------|
| cosm202 | 8882965 | z = 31 -> 0 in 890 steps, 6.9 h, one invocation; `[external noise] reading .../flamingo_panphasia_N2048/noise/` confirmed; `ZD_CornerModes = 1`, 2LPT; output 368 GB (312 GB checkpoint) |
| cosm203-206 | 8883893-8883896 | submitted 2026-10-01 |

## 6. Deferred: emulator-box runs (1215 Mpc/h, 4480^3)

- These can't share FLAMINGO's phases (section 1); they would use their own Panphasia descriptor.
- With the FLAMINGO-style descriptor (S1), `GridRes = 4480` gives relative level 12, N0 = 8192: 4097
  slabs of 1 GiB in double (about 4.4 TB per phase), generated once per phase and shared by all
  cosmologies. They would be cropped to 4480 (and converted to float, `EXTERNAL_NOISE_FLOAT`) before
  the runs.
- A descriptor with `S35` at relative level 6 would give N0 = 35 * 2^7 = 4480 exactly: no 8192 grid, no
  crop. Panphasia's support for S = 35 is unverified.
- Even S = 3 would give N0 = 6144 instead of 8192.

## 7. Validation ladder

| Test | What it checks | Where | Status |
|---|---|---|---|
| N = 256 tests 1-4 | loader, colouring, IC format (N0 = N) | `run256/validation_results.txt` | done (earlier) |
| Smoke test | production multistep reads external ICs, z = 31 to 0 | `smoke/` | 2LPT run (2026-09-30) reached the first LPT kick and aborted in `_pack_float`. These ICs have vel = displ exactly (f_cluster = 1), so `lpt_vel_scale` is double roundoff (7e-11), below the float32 kick's roundoff (1e-10). Production ICs (f_cluster < 1) have \|vel - displ\| ~ 3e-3 displ, like rockV2. Failed outputs are in `smoke/fail2_2lpt_pack/` and `...c202.fail2_2lpt_pack`. **S1b PASS** (`LagrangianPTOrder = 1`, `qsub -v LPT=1`, job 8881770): rc 0, z = 31 to 0 in 1583 steps (0.30 h on 1 node), slices at z = 30, 15, 7, 3, 1, 0.5, 0; external ICs unchanged (250/250 sha256 OK). |
| Round trip (a), N = 256 standalone | seed noise dump -> `LOAD_EXTERNAL_NOISE` reproduces the seed ICs, with qPLT, `Pk_norm = 8`, `NumZRanks = 2`; dumping does not perturb the RNG | `roundtrip256/results.txt` | **PASS** (2026-09-30): D and ICs bitwise identical (fraction 1.0, max ulp 0); plain vs seed+dump 0/501 files differ; <\|w\|^2> = 0.99981 |
| Crop (a), N = 300 from N0 = 512 | tests 1-4 plus test 1b (every mode of D is a positive real multiple of w_N0 at the same signed k) and r(k) with monofonIC's own crop | `crop300/run/validation_results.txt` | **PASS** (2026-09-30): test 1b phase max 5.9e-8 rad, fraction 1.0; r(k) = 1.00000 vs monofonIC's crop, zero translation; test 4 exact |
| Crop (b), N = 1800 from N0 = 2048 | all 1025 FLAMINGO slabs (<\|w\|^2> per plane, Hermitian y = 0), low-k identity with the N0 = 512 noise, per-plane k-space identity D = w L^-1.5 (flat P); the real-space block vs monofonIC's white noise is skipped (not written, see §3) | `crop1800/crop1800_results.txt` | **PASS** (2026-09-30): per-plane <\|w\|^2> in [0.965, 1.001], Hermitian 1e-15; low-k r = 1.00000 up to k = 128; D = w L^-1.5 to 6e-8 on 9 planes. The first check flagged 48 k-rows because `(fftfreq(N)*N).astype(int)` truncates 113.99999 to 113 at N = 1800 (0 cases at N = 256, 300). Fixed with `np.rint` in `crop_check_1800.py` and `validate_external_noise.py`. |
| Stage 5 | custom multistep with the real `FLAMINGO_nu_c202.par2` (only NP = 300^3, CPD = 125 and IC-step controls overridden) vs standalone `build_rt_load` from a hand-written par, same FLAMINGO noise: noise log line, `ZD_*` plumbing (sigma_ratio = 1, f_cluster), ICs to FFT roundoff (the builds plan FFTs differently, so not bytewise) | `stage5/results.txt` | Standalone half done (250 files). The first embedded try on 1 node passed the `ZD_*` plumbing check (22 keys, 0 differ), but `multistep` needs >= 2 X ranks with `NumZRanks = 2`. **PASS** on 2 nodes (job 8882857, 2026-10-01). The log reads the FLAMINGO noise; `ZD_*` keys equal the hand-written values; ICs agree to max\|diff\| 1.04e-6 rms (23 % bitwise, FFT plans differ). |
| Round trip (b), rockV2_n800 | embedded, production binary + seed vs custom binary + dumped noise; standalone seed vs embedded production | `roundtrip_rockv2/results.txt` | **PASS** (2026-09-30). The test run read the dumped noise; the ref run had no external-noise line; `ZD_*` match production except the (identical) `CLASS_power` path. Ref vs test: 74 % floats bitwise, max\|diff\| 9.5e-7 rms. Standalone seed vs embedded production: 21 % bitwise, 1.4e-6 rms (FFT plans differ). <\|w\|^2> = 0.99991. |

Round trip (a), Stage 5 and the crop tests gate the five production runs. If round trip (b) fails
later, stop runs still in progress and regenerate their ICs after the fix.

### Round-trip design (`DUMP_NOISE_SLABS`)

- A compile flag in `ZD_MPI_generation.c`. Where the seed path has set D (in double, before the float
  cast), it writes w = D / sqrt(P(k)) with the same `p->power(k)` that `zeldovich_ps_scaling` uses.
  The format is the loader's input format (`[x][z]` complex double per y-plane, N0 = N, masked modes 0).
- The seed path draws R = sqrt(-P log u), so <|w|^2> = 1, as for Panphasia.
- The hook only reads D, so the RNG stream is unchanged; the `plain` vs `seed` comparison checks this.
- Reconstruction gives sqrt(P) (D / sqrt(P)), equal to D to about 1 ulp in double. After the float cast
  the ICs are identical except for rare 1-float-ulp flips.
- FFT plans: `FFTW_MEASURE` can choose different plans run to run, so the standalone round trip builds
  all variants with `-DFFTW_PLANNER_FLAGS=FFTW_ESTIMATE`. D is compared first (before any FFT), then
  the IC files.

### rockV2_n800 round trip (b) (`roundtrip_rockv2/`)

- `pbs_seed_dump.sh` (`debug`, 1 node): standalone `build_rt_seed` with `param_rockv2_seed.par`, the
  exact ZD settings from the finished run's `abacus.par` (seed 12507, `Pk_norm = 8`, qPLT, 800^3,
  CPD 175, `NumZRanks = 2`). This writes the seed noise and a standalone copy of the ICs to
  `/flare/Abacus/helenshao/roundtrip_rockv2/seed_standalone/`.
- `pbs_embedded.sh` (`debug-scaling`, 4 nodes, 2 ranks per node, the same layout as production) runs the
  unchanged rockV2 par2 twice with `abacus.run -n 1 -P ICCleanupMode=Never -P
  InitialConditionsDirectory=<Flare> -P UseSCR=0`, once per binary:
  - `ref`: the production binary (eisenste 9575f04) with seed noise;
  - `test`: the custom binary, with `ZD_EXTERNAL_NOISE_DIR` pointing at the dumped noise.

  It then checks:
  - the `[external noise] reading` log line (present only in `test`);
  - the `ZD_*` keys against the production `abacus.par`;
  - `compare_rt256.py --fft-tol 1e-5` for `ref` vs `test`, plus the standalone seed ICs vs `ref`
    (standalone equals embedded at seed level).
- The two binaries plan their FFTs independently (`FFTW_MEASURE`), so the IC comparison allows FFT
  roundoff: max|diff| <= 1e-5 rms per field. The bitwise fraction is reported as well. When rank
  layouts differ (4 standalone vs 8 embedded ranks), particles are matched by (i, j, k).
