# External white noise (FLAMINGO / monofonIC-Panphasia) for massive-neutrino runs — cluster package

Branch: `neutrinos` (off `meson`). Goal: generate our Zel'dovich ICs from Willem Elbers' Panphasia
white noise (the FLAMINGO L1000N1800 realisation) coloured with **our own** CLASS \(P(k)\), and
validate it on a small N = 256 run before any production run.

Background, derivations and the Q&A behind every choice below: `NEUTRINOS_NOTES.md` (snapshot of
`Zeldovich-PNG/NEUTRINOS.md`; edit this copy on the cluster from now on).

## Contents

| File | What it is |
|---|---|
| `README.md` | This plan |
| `NEUTRINOS_NOTES.md` | Full notes: file format, half-cell phase, KK1025 fixing, `fixed_power`, Fix 4 mask, neutrino items, questions for Willem |
| `build_monofonic.sh` | Clone + build `wullm/monofonic@output_k_space_slabs` (no CLASS/zwindstroom/FastDF) |
| `monofonic_N256_test.conf` | monofonIC config: FLAMINGO descriptor and box at `GridRes = 256`, testing mode |
| `build_zmpi_extnoise.sh` | Meson build of this repo with `-DLOAD_EXTERNAL_NOISE=1 -DDUMP_D_SLABS=1` |
| `param_N256_extnoise.par` | Zeldovich-MPI parameters for the validation run |
| `run_N256_validation.sh` | Runs monofonIC → Zeldovich-MPI → validation, all output to log files |
| `validate_external_noise.py` | Tests 1–3, numbers written to `validation_results.txt` |

## What the code does (branch `neutrinos`)

- `config.h`: `LOAD_EXTERNAL_NOISE`, `EXTERNAL_NOISE_DIR` (`"external_noise"`), `EXTERNAL_NOISE_PREFIX`
  (`"output_k_space_slab_y."`), `EXTERNAL_NOISE_FLOAT` (0 = double), `DUMP_D_SLABS`, `DUMP_D_DIR` (`"D_slabs"`).
  Paths are relative to the directory the executable is run from.
- `load_external_noise_slab()` (`ZD_MPI_generation.c`): reads slab \(k_y\) = `global_y`, infers \(N_0\) from the
  file size (even, \(\ge N\)), crops by signed wavenumber when \(N_0>N\), transposes `[x][z]` → our `[z][x]`,
  `MPI_Abort` on a missing/short file. Each rank reads only the slabs it generates.
- In the z-loop, modes that pass the existing mask (not DC, not Nyquist, inside the cutoff sphere if
  `CornerModes = 0`) get \(D = \sqrt{P(k)}\,w\) via `zeldovich_ps_scaling()`; no RNG draw.
  `fixed_power` is ignored (the file already has FLAMINGO's KK1025 partial fixing).
- `DUMP_D_SLABS`: writes the input \(D\) of each primary slice to `D_slabs/D_slab_y.<y>`
  (`[x][z]`, complex double, same layout as Willem's files) for tests 2 and 3.
- Everything is `#if`-guarded; with the flags at 0 the code is the `meson` Gaussian generator.

## Step-by-step on the cluster

### 0. Get the code
```bash
cd <Zeldovich-MPI clone>
git fetch origin && git checkout neutrinos && git pull
python3 -c "import numpy, h5py"            # needed by the validation script
```

### 1. Build monofonIC (once)
```bash
module load <cmake mpi fftw gsl hdf5 fortran modules>
MONO_DIR=$HOME/monofonic bash neutrino_validation/build_monofonic.sh > build_monofonic.log 2>&1
```
Needs FFTW3 double with MPI + OpenMP, GSL, HDF5, a Fortran compiler. Register for the Panphasia licence
(<http://icc.dur.ac.uk/Panphasia.php>) before using these ICs in anything published.

### 2. Build Zeldovich-MPI with the external-noise flags
```bash
module load frameworks fftw/3.3.10
bash neutrino_validation/build_zmpi_extnoise.sh > build_zmpi_extnoise.log 2>&1
# -> build_extnoise/src/Zeldovich_MPI
```

### 3. Prepare the P(k) file
- Two columns: \(k\) [h/Mpc], \(P(k)\) [(Mpc/h)\(^3\)], `#` comments allowed.
- At `InitialRedshift` (31 in the template) and already the spectrum you want applied
  (for neutrino runs: \(P_{cb}\), see notes §4). `ZD_Pk_norm = 0` means no \(\sigma_8\) rescaling:
  the code uses the file as is (`p->power(k) = P_file / BoxSize^3`).
- Extend to \(k \ge k_{\rm Nyq}\sqrt{3}\) = 2.07 h/Mpc for N = 256, L = 673.2 (and ≈ 14.5 h/Mpc for
  N = 1800), otherwise the spline extrapolates (the code warns).

### 4. Run the validation (interactive compute node)
```bash
MONO_EXE=$HOME/monofonic/build/monofonIC \
ZMPI_EXE=$PWD/build_extnoise/src/Zeldovich_MPI \
PK_FILE=/path/to/pk_cb_z31.txt \
bash neutrino_validation/run_N256_validation.sh
```
Options: `RUN_DIR` (default `neutrino_validation/run256`), `MONO_RANKS`, `ZMPI_RANKS`, `NTHREADS`,
`MPIEXEC`, `MPI_OPTS` (default is Aurora-style `--ppn --depth --cpu-bind depth`), and
`SKIP_MONO=1` / `SKIP_ZMPI=1` to redo only later steps.

Outputs in `run256/`:

| Path | From |
|---|---|
| `mono/monofonic.log`, `mono/mono_N256_delta.hdf5` | step 1 |
| `external_noise/output_k_space_slab_y.{0..128}` | step 1 (129 files × 1 048 576 bytes) |
| `zmpi_run.log`, `zmpi_out/`, `D_slabs/D_slab_y.{0..128}` | step 2 |
| `validate.log`, `validation_results.txt` | step 3 |

### 5. Pass criteria

| Test | Pass |
|---|---|
| 1. White noise (slabs only) | \(\langle\|w\|^2\rangle \approx 1\) in every shell; \(\|w\| = 1\) exactly for \(k^2 \le 1025\); Nyquist = 0 |
| 2. \(P_{\rm ours}/P_{\rm input}\) (D slabs + P file) | ratio ≈ 1 in all shells, scatter shrinking with \(k\) |
| 3. vs monofonIC δ (D slabs + HDF5) | \(r(k) \approx 1\); best-fit translation ≈ (0,0,0) cells; small residual phase |

### 6. If something fails

- **Test 2 off by a constant.** \(N^6\) or \((2\pi)^3\): Fourier convention; try `--ours-norm` (V/N⁶).
  \(h^3\): units of the P file. \(\sigma\)-ratio²: `ZD_Pk_norm` was not 0.
- **Test 2 wrong shape.** \(k\) units (BoxSize in Mpc/h ↔ \(k\) in h/Mpc; `ZD_Pk_scale`).
- **Test 3 \(r\) low everywhere.** Axis order: our transpose, or monofonIC's HDF5 order. Transpose the HDF5
  array in the script (`d = d.transpose(...)`) before suspecting the reader.
- **Test 3 \(r \approx 1\) but shift ±0.5 cell.** Half-cell phase convention (notes §2.1); decide whether to correct.
- **Test 3 \(r\) drops only near \(k_{\rm Nyq}\).** Expected: Nyquist/corner masking and monofonIC's own filtering.
- **`ERROR: ... not an even N0 x N0 grid`.** Wrong precision (`EXTERNAL_NOISE_FLOAT`) or truncated file.

### 7. After the validation passes

1. Abacus run from the N = 256 ICs (`zmpi_out/`) as a smoke test.
2. Production, FLAMINGO-matched: monofonIC `GridRes = 1800` → \(N_0 = 2048\), ≈ 69 GB of slabs;
   our N = 1800 crops to \(|k_i| \le 900\) = FLAMINGO's modes. Build with `DUMP_D=0`
   (the dump would be ≈ 47 GB). Or run monofonIC once and keep the slabs on the project filesystem.
3. Settle the neutrino items (notes §4): \(P_{cb}\) vs \(P_m\), gauge, scale-dependent \(f(k)\) for velocities,
   LPT order (monofonIC FLAMINGO uses 3LPT), `ZD_f_cluster`, neutrino particles, units.
4. Questions for Willem: notes §6 (half-cell phase, `mode_weightings` placement, float vs double files,
   whether the N = 1800 FLAMINGO ICs used exactly this descriptor).

## Reconciling with branch `PNG` later

`PNG` (commit `c9263ae`, local only) holds Grace's PNG blocks plus an earlier copy of this mode. When Grace's
missing pieces arrive: merge `neutrinos` into `PNG`. The expected conflict is the in-loop block — on `PNG`
the colouring needs its own mask check because Grace's load blocks overwrite \(D\) after the mask.
