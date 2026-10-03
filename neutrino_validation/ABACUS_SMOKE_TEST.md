# Abacus smoke test for the external-noise ICs (N = 256, cosm202)

Goal: show that the Zeldovich-MPI ICs built from Willem's Panphasia white noise (`run256/zmpi_out/`)
are read by Abacus correctly and then evolve as expected. The checks go from exact identities
(bookkeeping, IC round-trip) to linear-theory growth, and only then to non-linear statistics.
Each stage has pass criteria, writes its numbers to a text file, and gates the next stage.

Not in scope: the neutrino physics decisions (P_cb vs P_m at z_start, gauge, scale-dependent f(k)
for velocities, LPT order; `NEUTRINOS_NOTES.md` §4) and production resolution. The smoke test
checks that the pipeline is correct for the physics we chose, not that the physics choice is right.

## 0. Decisions and prerequisites

### 0.1 Do we need to recompile anything?

| Path | Zeldovich-MPI | Abacus |
|---|---|---|
| **A. External ICs** (stages 1–4): Abacus reads `zmpi_out/` with `ExternalICs = 1` | No. The standalone `build_extnoise/src/Zeldovich_MPI` already made the ICs | No rebuild against `neutrinos` is needed: with `ExternalICs = 1` the embedded generator never runs (`multistep.cpp`, `P.EmbeddedZeldovichMPI()`). We do need a working `multistep`, though, and the one we have (`abacus/build/src/Multistep/multistep`, Aug 9) is a **debug** build of an older checkout. Rebuild in release mode on a compute node |
| **B. Embedded ICs** (stage 5): `multistep` generates the ICs itself from the noise | Yes: `abacus/subprojects/zeldovich_mpi` must be on `neutrinos` | Yes: reconfigure with `-Dzeldovich_mpi:extra_cpp_args=-DLOAD_EXTERNAL_NOISE=1` (no `DUMP_D_SLABS`) and rebuild |

**Decision (2026-09-27):** path A uses the production build instead of a rebuild:
`/home/eisenste/abacus-store/9575f04bdaf2f7c88c2b160ac0d60b428a80b4d3` (env `env/aurora.sh`, `-O3`,
built 2026-09-14), the same code as the AbacusAurora emulator runs. Derivatives for CPD = 125
(`deriv32_125_8_2_8`) already exist in `/flare/Abacus/Derivatives`.

Submodule provenance: the parent repo records `zeldovich_mpi` at `c46448a` (`meson`), but the working tree
is at `e3bf1a7` (`neutrinos`). With the flags off, `neutrinos` builds the same generator as `meson`
(every new block is `#if`-guarded), so a path-A rebuild with `neutrinos` checked out is harmless. Record
`git -C subprojects/zeldovich_mpi rev-parse HEAD` in the run provenance either way.

Rebuild (compute node, see `.cursor/rules/aurora-meson-fftw.mdc`):
```bash
qsub -I -l select=1 -l walltime=1:00:00 -q debug -A Abacus -l filesystems=home:flare
cd /home/helenshao/InitialConditions/abacus && source env.sh
meson setup build_smoke --buildtype=release        # path B: add -Dzeldovich_mpi:extra_cpp_args=-DLOAD_EXTERNAL_NOISE=1
meson compile -C build_smoke > build_smoke.log 2>&1
```

### 0.2 Fixed facts the configuration must respect

| Item | Value | Why |
|---|---|---|
| `NP`, `BoxSize` | `256**3`, `673.2` | ICs |
| `CPD` | **125** | `zmpi_out` has one `ic_NNNN` per x-slab for `CPD = 125` (125 files). Another CPD means regenerating the ICs (cheap) and rerunning test 4 |
| `NumZRanks` | 1 | flat `ic_*` layout (Zeldovich mode 3) |
| `ICFormat` | `"RVZel_2D"` (same as `"RVZel"` when `NumZRanks = 1`) | 32-byte record: 3 × `uint16` (i, j, k), padding, `float32` displ[3], vel[3] |
| `ICPositionRange` | `673.2` (= BoxSize) | test 4: `-i k·psi = +1.0000000 × FFT(dens)`, so displacements are in Mpc/h |
| `ICVelocity2Displacement` | `1.0` | test 4: `vel = 1.0000000 × displ` exactly (`f_cluster = 1`) |
| `InitialRedshift` | `31.0` | P(k) file is P_cb(z = 31), used as-is (`ZD_Pk_norm = 0`) |
| Axis map | Abacus (x, y, z) = lattice (i, j, k) = D-slab / monofonIC axes (**z, y, x**) | test 4 (b). Swapping two axes is a mirror image, so statistics are unchanged, but any mode-by-mode or halo-by-halo comparison with monofonIC/FLAMINGO must swap x ↔ z |
| Cosmology | `abacus_cosm202/cosm.def` (FLAMINGO Planck, `Omega_Smooth = 0.001422`) | Abacus treats the 60 meV neutrinos as a smooth, non-clustering component |

### 0.3 Things that would silently break the test

- **IC deletion.** After the IC step, `CleanupICs()` with `MPI_size_z == 1` makes rank 0 run `remove_all`
  on the whole `InitialConditionsDirectory`. With `ExternalICs = 1` the default `"Auto"` would not delete
  anything, but set `ICCleanupMode = "Never"` anyway, and point Abacus at a **copy** of `zmpi_out`,
  never at `run256/zmpi_out` itself.
- **`base.def` defaults that don't apply here:** `InitialRedshift = 99`, `ZD_Pk_norm = 8.0`, `ZD_qPLT = 1`,
  lightcones, MapLog, MAPs. Override all of them in the smoke `.par2`.

## Stage 0 — IC gate (done, 2026-09-27)

`validate_external_noise.py` tests 1–4 on `run256/`. Outputs: `validation_results.txt` (tests 1–3),
`validation_test2_per_mode.txt`, `validation_test4.txt`.

| Test | Result |
|---|---|
| 1 white noise | `<|w|^2>` = 0.992 overall, 0.994–1.013 per shell; fixed modes \|w\| = 1 to 3e-16; Nyquist = 0 |
| 2 colouring (per mode) | `P_ours / (P_in |w|^2)` − 1 ≤ 7e-5; phase(D) = phase(w) to 6e-8 rad |
| 3 vs monofonIC (`--mono-sign -1`) | r(k) = 0.99997–1.00000; translation (0, 0, 0); residual phase 1.2e-8 rad |
| 4 particle files | 0 missing / 0 duplicate sites; displ[c] ∥ axis c (off-diagonal ≤ 2e-4); curl/div = 5e-7; FFT(dens) = N³ D and −i k·psi = N³ D, rms per-mode residual ≤ 2.6e-6 (single-precision FFT pipeline; the float32 storage floor alone is ~1e-7) |

## Stage 1 — Staging (login node is fine; plain file copy)

```bash
SMOKE=/flare/Abacus/helenshao/smoke_nu256_c202          # or another flare path
mkdir -p $SMOKE/ic && cp -p run256/zmpi_out/ic_* run256/zmpi_out/dens_* $SMOKE/ic/
(cd run256/zmpi_out && sha256sum ic_* dens_*) > $SMOKE/ic_sha256_source.txt
(cd $SMOKE/ic && sha256sum -c ../ic_sha256_source.txt) > $SMOKE/ic_sha256_check.txt
chmod a-w $SMOKE/ic/*                                     # Abacus only reads them
```
Pass: every line of `ic_sha256_check.txt` says `OK`; 125 `ic_` + 125 `dens_` files, 536 870 912 bytes of `ic_*`.

**Done (2026-09-27)**, at `/flare/Abacus/helenshao/ICs/nu256_c202_panphasia/ic` instead: outside every
directory `abacus.run --clean` wipes (Output/Working/State), with the directory itself also `a-w`.
125 + 125 files, 536 870 912 bytes, all checksums `OK`.

## Stage 2 — Parameter file `smoke_nu256_c202.par2`

```text
#include "./$ABACUS_SITE$.def"
#include "../Cosmologies/abacus_cosm202/cosm.def"
#include "./base.def"

SimName = AbacusAurora_smoke_nu256_c202
BoxSize = 673.2
NP = 256**3
CPD = 125
NumZRanks = 1
InitialRedshift = 31.0

ExternalICs = 1
UseLegacyZeldovich = 0
InitialConditionsDirectory = "/flare/Abacus/helenshao/smoke_nu256_c202/ic"
ICCleanupMode = "Never"
ICFormat = "RVZel_2D"
ICPositionRange = @BoxSize@
ICVelocity2Displacement = 1.0
LagrangianPTOrder = 2          # production setting; run LagrangianPTOrder = 1 as variant S1b

# Full particle outputs: an early one for the IC round-trip, then linear-regime and late times
TimeSliceRedshifts = [30.0, 15.0, 7.0, 3.0, 1.0, 0.5, 0.0]
OutputFormat = "RVdoublePID"   # default RVdouble has no PIDs; stage 4 needs (i, j, k) per particle
FinalRedshift = 0.0

# Off for the smoke test
LCOrigins = None
LCIsOctant = None
LCZMax1 = None
LCZMax2 = None
LCBoxRepeats = 0
LCSamplerInner = None
LCSamplerOuter = None
LCDoHealpix = None
LCHealpixZMax = None
MapLogRedshifts = None
TimeSliceRedshifts_Subsample = None
TimeSliceRedshifts_Samplers = None
LogVerbosity = 0
```
Check the parsed result before submitting: `python -m abacus.param smoke_nu256_c202.par2 -o smoke_parsed.par > param_parse.log 2>&1`,
then grep `smoke_parsed.par` for every key in §0.2. Each `None` override has to be accepted by the parser;
if one isn't, delete that key instead. PPC = 256³/125³ ≈ 8.6 is low but fine here.

**Done (2026-09-27):** `AbacusAurora/Simulations/smoke_nu256_c202.par2` (the draft above plus `FlipZelDisp = 0`,
`UseSCR = 0`, `StateIOMode = "overwrite"`, `CompressWriteState = 0`, `SCR_Checkpoint_Seconds = 0`,
`ExitAfterCheckpoint = 0`, and the staged IC path). Parsed with the production env into `smoke/smoke_parsed.par`:
all `None` keys accepted, every §0.2 key as intended, `mpirun_cmd` = 2 ranks on 1 node, outputs in
`/flare/Abacus/helenshao/AbacusAurora_smoke_nu256_c202`. Leftover `ZD_*` keys are ignored with `ExternalICs = 1`
(`preprocess_params` and `make_ic` both skip).

**PID and output format (confirmed in the code):** the RVZel loader calls `aux.init_dm_aux(ijk)` → `set_pid(i, j, k)`,
and `pid()` = `i | j << 16 | k << 32`. `RVdoublePID` writes 56 bytes per particle: `double pos[3]` in unit-box
units centred on 0 (the lattice site is `ijk/ppd − 0.5`, `ZelPos` in `lpt.cpp`), `double vel[3]` in unit-box
redshift-space units, `uint64 tag` = `pid()`. `FlipZelDisp` (default 0) would negate the displacements.

## Stage 3 — Run and live checks (one debug node)

Script: `smoke/pbs_smoke_run.sh` (`qsub` it; logs `pbs_smoke_run.log`, `run.log`, `provenance.txt`,
`ic_sha256_pre.txt`/`ic_sha256_post.txt` in `smoke/`). It refuses to start if the output directory already exists.

`python -u -m abacus.run smoke_nu256_c202.par2 > run.log 2>&1` inside a PBS job (1 node, `debug`,
1 h; if the MPI build insists on more than one rank, use 2 nodes with `NumZRanks = 1`). Everything goes to `run.log` plus Abacus's own logs.

| Check | Where | Pass |
|---|---|---|
| Embedded generator skipped | log: `Embedded Zeldovich_MPI disabled (ExternalICs=1 ...) loading pre-existing ICs from ...` | present |
| All particles read | IC step summary / `NP` in state | 16 777 216 |
| IC files untouched | re-run the stage-1 `sha256sum -c` after the job | all `OK` |
| 2LPT consistency | `lpt_vel_scale` (max \|vel·Canonical_to_VelZSpace − displ\|) in the IC/LPT log | ≈ 0 up to unit rounding, since vel = displ |
| Stability | every step: no NaN/Inf, `max accel` and time step smooth, no QUIT/assert | reaches z = 0 |
| Resources | step time, memory high-water mark | recorded for scaling to N = 1800 |

## Stage 4 — Quantitative acceptance (compute node; outputs in `$SMOKE/analysis/*.txt`)

Every check works mode by mode, using the fixed phases and, below k² ≤ 1025, the fixed amplitudes.
That makes it much sharper than a shell-averaged P(k) comparison. Particle IDs give each particle's (i, j, k),
so displacements can be rebuilt exactly on the lattice (same machinery as test 4). Before relying on this,
confirm how the `RVdoublePID` PID encodes (i, j, k); `aux.init_dm_aux(ijk, ...)` in `loadIC.cpp` sets it.

**S4.1 IC round-trip (z = 30 slice).** Rebuild psi(q) = x − q (with periodic wrapping) and compare it per mode
with the IC displacement:
- correlation r > 0.99999 for k < k_Nyq/2;
- amplitude ratio equal to the 2LPT-corrected linear growth D(30)/D(31), within 1e-3;
- no axis permutation or sign flip (the test-4 matrix computed on the Abacus output has to stay diagonal and +1).

This checks that Abacus read positions, axes, units and velocities the way test 4 assumed.

**S4.2 Linear growth per mode (z = 15, 7, 3).** For each k < 0.05 h/Mpc, fit g(k, z) = δ_sim(k, z) / δ_IC(k). Compare it with:
- (a) Abacus's own growth table for cosm202 (with `Omega_Smooth`), from its log or from `abacus.cosmo`: agree to ≤ 2e-3 at k ≤ 0.03 h/Mpc;
- (b) CLASS, sqrt(P_cb(k, z) / P_cb(k, 31)). Expect Abacus to track CLASS's small-scale limit: CLASS's own P_cb growth
  from z = 31 to 0 is 24.880 for k ≥ 0.3 h/Mpc, rising to 24.980 at the fundamental mode k = 0.0093 h/Mpc.
  That 0.4 % scale dependence comes from neutrino clustering plus synchronous-gauge near-horizon terms, which smooth
  neutrinos and Newtonian gravity don't reproduce. Record g_sim / g_CLASS(k); a trend of that size is expected and
  documents the gauge/neutrino item in `NEUTRINOS_NOTES.md` §4, not a bug.

**S4.3 Bookkeeping at every slice.** Particle count = 256³; IDs unique and covering the lattice; momentum drift ≈ 0.

**S4.4 Non-linear sanity (z = 1, 0.5, 0).** P(k) against CLASS linear × growth on large scales and against an emulator/halofit shape;
halo mass function if group finding is turned on in a follow-up run. Compare with a **Gaussian control**: the same
Zeldovich-MPI build without `LOAD_EXTERNAL_NOISE` (`meson` generator), same P(k) file, same par2 apart from the IC directory.
The two runs must agree within cosmic variance; the fixed-noise run should scatter less at k² ≤ 1025.

**S4.5 Optional: FLAMINGO large scales.** If a FLAMINGO L1000N1800 DMO z = 0 density field is available, cross-correlate
at k < 0.1 h/Mpc after swapping x ↔ z (§0.2) and converting Mpc → Mpc/h (L = 1000 Mpc = 673.2 Mpc/h, h = 0.6732 vs our 0.673).
Expect r(k) → 1 on the largest scales. Different resolution, LPT order (FLAMINGO 3LPT) and neutrino treatment
limit the agreement at higher k. Before relying on the swap, confirm FLAMINGO's particle axis order with Willem.

Go/no-go: S4.1–S4.3 must pass exactly as stated. S4.2(b) and S4.4 are recorded and judged. A failure in S4.1
means stopping and re-checking the §0.2 conventions before looking at anything later.

## Stage 5 — Embedded path (production mechanism)

1. Rebuild Abacus as in §0.1 path B.
2. Stage the slabs where `multistep` looks for them: `EXTERNAL_NOISE_DIR` (`"external_noise"`) is relative to
   the **working directory of the multistep process**. Confirm which directory `abacus.run` launches it in
   and symlink `external_noise/` there.
3. Override in the par2: `ExternalICs = 0`, `UseLegacyZeldovich = 0`, `ZD_Pk_filename` = the z = 31 P_cb file,
   `ZD_Pk_norm = 0` (`base.def` has 8.0), `ZD_qPLT = 0` (to match `zmpi_out`), `ZD_CornerModes = 1`, `ZD_k_cutoff = 1`,
   `ZD_qPk_fix_to_mean = 0`, `ZD_f_cluster = 1.0`, `ICCleanupMode = "Never"`. Also check whether the embedded path
   applies any growth rescaling from `ZD_Pk_file_redshift` (in `cosm.def`) when `ZD_Pk_norm = 0`.
4. Stop after the IC step (e.g. `FinalRedshift` just below 31, or a HALT file) and compare the generated `ic_*`
   with `zmpi_out`: they should be byte-identical (same build flags and precision). If they aren't, run test 4
   on them (`--ic-dir`) and require the same numbers as stage 0.

## Stage 6 — Toward N = 1800

- Cropping (N0 > N) is still untested. Run monofonIC with `GridRes = 300` (N0 = 512) and Zeldovich-MPI with N = 256,
  then tests 1–4. Test 3 must still give r = 1 and zero translation after the crop.
- Test 4 as written loads everything into one process (4.1 GB RSS, 55 s at N = 256). At N = 1800 it needs a
  sub-volume or MPI version, or per-slab checks of the identities that are local in x.
- Then decide the N = 1800 CPD and redo stages 1–4 at reduced output.

## Open items carried by this test

- P_cb at z = 31 taken straight from CLASS (synchronous gauge) vs back-scaling from z = 0 (monofonIC/zwindstroom).
- A single f for velocities (vel = displ) vs scale-dependent f(k) with massive neutrinos.
- 1LPT(+2LPT in Abacus) at z = 31 vs FLAMINGO's 3LPT.
- Single-precision IC pipeline (per-mode residual ~1e-6); the `double-precision` branch could lower it if needed.
