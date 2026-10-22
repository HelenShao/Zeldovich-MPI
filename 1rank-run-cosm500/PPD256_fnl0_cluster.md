# PPD = 256, f_NL = 0 sanity check

Compare ordinary Zel'dovich initial conditions (meson branch, no PNG generation) with the PNG load pass at `ZD_f_NL = 0`. The particles should match to roundoff. This run is 1 MPI rank, `NP = 256^3`.

Use `project1/Zeldovich-PNG-grace` for the PNG side. The `png-port` branch of `Zeldovich-MPI` only has parameter getters. It does not write `Phi_gaussian_y.*` or read `kernel_y.*`.

Pull `my-edits-on-meson` on the cluster. That branch has Grace's C kernel port (`src/kernel`, commit `e0ad77b`) and the y-keyed slice names (`ZD_PNG_DIR/Phi_gaussian_y.{y}`, `kernel_y.{y}`). This note uses the Python kernel and the Zel'dovich dump/load path. Grace's `Kernel_MPI` binary reads and writes the same y-keyed slices, so it can replace the Python kernel in step 2.

Build on a compute node (or a compile node) that has meson, ninja, MPI, single-precision FFTW with OpenMP, and bison 3 or newer. Keep `png_slices/` and both initial-condition directories on a filesystem shared with the node that runs `kernel_MPI.py`.

Save every log named below. `run.sh` already writes its own stdout log. For the Python kernel and the final position/velocity comparison, redirect stdout and stderr into a text file yourself.

## One-time edits in this bundle

Edit `cosm500_abacus_legacy_zeldovich.par2`:

- `NP = 16777216` (256³). It is currently `2097152` (128³).
- `CPD = 256`. It is currently `243`, and 243 does not divide 256. The writer uses `file_index = z * CPD / N` and puts `N/CPD` planes in each file.
- `ZD_f_NL = 0`. The dump pass ignores this value. Setting it here keeps the later load pass from applying `f_NL = 100`.
- `NumZRanks = 1`. The template is currently `2` (Grace's 128³, 2-rank cosm500 run). This sanity check is 1 rank. `ZD_NumZRanks` is no longer in the file.

Edit `run.sh`: the wisdom preflight is `N=128`. Change that assignment to `N=256`. The particle grid still comes from `NP`; this `N` is only the size passed to `wisdom_rank0`.

`build.sh` and `run.sh` default `REPO_ROOT` to two directories above this bundle (`project1/`). Point it at the directory that contains `meson.build`:

```bash
export REPO_ROOT=/path/to/Zeldovich-PNG-grace
cd "$REPO_ROOT/1rank-run-cosm500"
```

Leave `BoxSize = 325` and the rest of the cosmology, seed (`ZD_Seed = 12345`), `ZD_Pk_norm = 8`, `ZD_Pk_sigma`, `ZD_k_cutoff = 1`, and `ZD_qPLT = 1` as they are. Both codes must see those same values.

`data/class_pk_z0.dat`, `data/class_pk_primordial_dimensional.dat`, and `data/eigmodes128` are already in this bundle. `env.sh` picks them up if you do not override `PK_FILE` and `PLT_FILE`.

## Shared settings for the meson run

Check out the `meson` branch of `Zeldovich-MPI` (no PNG generation). Give it a parameter file with the same box, `NP`, `CPD`, seed, `ZD_Pk_filename` (`class_pk_z0.dat`), PLT file, `ZD_Pk_norm`, `ZD_Pk_sigma`, `ZD_Pk_smooth`, `ZD_k_cutoff`, `ZD_qPLT`, `ZD_Version`, `InitialRedshift`, `ICFormat = RVZel`, and `NumZRanks = 1`.

Omit these lines from the meson file. That branch does not install them, and the parser will abort on an unknown name:

- `ZD_PNG_DIR`
- `ZD_Pk_primordial_filename`
- `ZD_f_NL`, `Omega_M`, `ZD_n_s`
- `ZD_f_NL_local`, `ZD_f_NL_equil`, `ZD_f_NL_orthoCMB`, `ZD_f_NL_orthoLSS`

Write the meson particles somewhere that will not be overwritten by the PNG runs, for example `IC_OUTPUT_DIR=.../output_meson_ppd256`.

## 1. Dump pass (write `Phi_gaussian_y.*`)

`config.h` defaults to `DUMP_GAUSSIAN_D=1` and `LOAD_D_FROM_FILE=0`. Leave `build.sh` at `-Dextra_cpp_args=-g` for this pass.

```bash
export REPO_ROOT=/path/to/Zeldovich-PNG-grace
export MESON_BUILD="$REPO_ROOT/build/ppd256_dump"
export IC_OUTPUT_DIR=/path/to/output_png_dump_ppd256
export PNG_DIR=/path/to/png_slices
cd "$REPO_ROOT/1rank-run-cosm500"
./run.sh
```

Stdout goes to `cosm500_abacus_legacy_zeldovich_bin_generation.log` in this bundle. This binary writes `png_slices/Phi_gaussian_y.{y}` for `y = 0 .. 128` and also writes Gaussian particle initial conditions under `IC_OUTPUT_DIR`. It does not open kernel files.

Do not delete `png_slices/` after this step.

## 2. Build the kernels

Run this from `1rank-run-cosm500` so the default `PARAM_FILE=.runtime/param.par2` resolves. That file is what `run.sh` wrote, with absolute paths. Use `kernel_MPI.py`. `kernel_MPI_test.py` still has the old rank-keyed names.

```bash
cd "$REPO_ROOT/1rank-run-cosm500"
python3 kernel_MPI.py > kernel_ppd256.log 2>&1
```

The script reads `Phi_gaussian_y.0` to get `N`, then writes `png_slices/kernel_y.{y}`. The kernel grid is `512^3`, so give the process a few GB of memory. Confirm `kernel_ppd256.log` prints `N 256` and that `kernel_y.0` through `kernel_y.128` exist.

To use the C kernel instead, build it in double precision and run it on the same parameter file. It reads `ZD_PNG_DIR` from the par2 file and writes `kernel_y.{y}` into that directory. Any rank count works, as long as it is a multiple of `NumZRanks` and at most `N/2 + 1`. It does not have to match the dump or load runs.

```bash
cd "$REPO_ROOT"
tests/kernel_regression/build_kernel.sh
cd "$REPO_ROOT/1rank-run-cosm500"
mpiexec -n 2 "$REPO_ROOT/build/kernel_double/src/Kernel_MPI" .runtime/param.par2 > kernel_ppd256.log 2>&1
```

Check the same things: `N=256` in the log and `kernel_y.0` through `kernel_y.128` in `png_slices/`.

## 3. Load pass (`f_NL = 0`)

This is a different binary. Passing only `-DLOAD_D_FROM_FILE=1` leaves the header default `DUMP_GAUSSIAN_D=1` in place, and the mutual-exclusion `#error` fires. In `build.sh`, change the meson line to:

```bash
-Dextra_cpp_args="-DDUMP_GAUSSIAN_D=0,-DLOAD_D_FROM_FILE=1,-g"
```

Use a second build directory and a second initial-condition directory. Keep the same `PNG_DIR`.

```bash
export MESON_BUILD="$REPO_ROOT/build/ppd256_load"
export IC_OUTPUT_DIR=/path/to/output_png_load_fnl0_ppd256
export PNG_DIR=/path/to/png_slices
cd "$REPO_ROOT/1rank-run-cosm500"
./run.sh
```

`ZD_f_NL` is read from the parameter file at runtime. With the template set to 0, this run reads `kernel_y.{y}`, adds `0 * K`, and the `3/5` and `5/3` factors cancel. The particles under `output_png_load_fnl0_ppd256` are the ones to compare to `output_meson_ppd256`.

Put the `build.sh` extra-args line back to `-Dextra_cpp_args=-g` when you want another dump binary. The two `MESON_BUILD` directories keep the dump and load executables from overwriting each other.

## 4. Compare

Compare positions and velocities in `output_meson_ppd256` with `output_png_load_fnl0_ppd256` (same `RVZel` reader on both). Write the maximum absolute position difference and the maximum absolute velocity difference to a text file, for example `fnl0_ic_diff.txt`. A match at roundoff means the PNG round trip at `f_NL = 0` reproduces the meson Zel'dovich field.

The dump directory `output_png_dump_ppd256` is a second Gaussian copy from the PNG tree. A mismatch that appears against meson but not against the dump directory is a driver or parameter difference. A mismatch against both is in the load path.

Full time evolution with `abacus.py` is the bispectrum test, after this initial-condition comparison has passed. Point it at `InitialConditionsDirectory` only once `fnl0_ic_diff.txt` is at roundoff. Then raise `ZD_f_NL` (the template value used to be 100) and rerun only the load binary. The kernel files stay the same.
