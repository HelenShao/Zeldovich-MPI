# Kernel_MPI regression tests

`Kernel_MPI` (`src/kernel/`) is the MPI/C port of `1rank-run-cosm500/kernel_MPI.py`.
These tests check it against that confirmed Python code and against an independent
plain-numpy implementation.

## Build and run

```bash
tests/kernel_regression/build_kernel.sh                 # -> build/kernel_double/src/Kernel_MPI (double precision)
python3 tests/kernel_regression/run_regression.py       # quick tier, about 1 minute
meson test -C build/kernel_double kernel_regression_quick   # same quick tier, via meson
python3 tests/kernel_regression/run_regression.py --tier full   # real dumps in 1rank-run-cosm500/bin_files
```

The Python interpreter must have numpy and pyfftw. The default is the interpreter running the script; `--python` overrides it.
Useful options: `--fnl equil` and `--N 16` (restrict cases, repeatable), `--omp-threads`, `--oversubscribe`,
`--json results.json`, and `--work-dir`. The default work dir is `build/kernel_regression_work`, where every
run's outputs and logs are kept. Exit status is 0 if all checks pass, 1 on any failure, and 77 if the run was skipped.

## What is checked

| tag | comparison | gated |
|---|---|---|
| `py` | C (compat mode, M=2N) vs `kernel_MPI.py`, same rank count | yes |
| `direct` | C vs `reference_direct.py` (full-field numpy, no rank emulation) | yes |
| `ranks` | C on each rank layout vs C on 1 rank (matched by y, not by file); forced staged Alltoallv vs normal | yes |
| `nyq` | C with `KERNEL_ZERO_NYQUIST=1` vs `reference_direct.py --zero-nyquist` | yes |
| `pad` | zero-Nyquist C with M=3N/2 vs M=2N (both alias-free) | yes |
| `nyq`, `pad` INFO | size of the Nyquist / packing effect in compat mode | no |

Gated tolerances: relative L2 ≤ 1e-10 and max relative error ≤ 1e-8 over modes above 1e-8·max.
The y=0 output plane must also be Hermitian away from the Nyquist row and column.
The quick tier uses synthetic Hermitian inputs (`make_test_phi.py`) with N ∈ {16, 32} and rank layouts
(R, NumZRanks) ∈ {(1,1), (2,1), (2,2), (3,1), (4,2), (6,2)}, including uneven pair splits.
The Python reference outputs are cached by a hash of the inputs and of `kernel_MPI.py`.

To compare any two output directories by hand:

```bash
python3 tests/kernel_regression/compare_kernel_outputs.py <dir_under_test> <reference_dir>
```

## Kernel_MPI runtime switches

| variable | effect |
|---|---|
| `KERNEL_PYTHON_COMPAT=1` (default) | default M = 2N, as in `kernel_MPI.py` |
| `KERNEL_PYTHON_COMPAT=0` | default M = smallest valid M ≥ 3N/2 |
| `KERNEL_PAD_M=<M>` / par2 `KernelPadM` | explicit padded size (M % 4 == 0, M − N even) |
| `KERNEL_ZERO_NYQUIST=1` | zero the kx/ky/kz = N/2 modes of Φ before the kernels |
| `KERNEL_ONLY_TERM=<A1…G1>` | debug: run one kernel shape with coef 1 (output = J·Q) |
| `KERNEL_FORCE_STAGED_A2A=1` | test: always use the INT_MAX-safe staged Alltoallv path |

## Known difference in compat mode: Nyquist modes and the packing trick

`kernel_MPI.py` (with `PACK_GH=True`) and compat mode evaluate g·h as Im(z²)/2, where z = IFFT[(G + iH)Φ].
This is exact only if g and h are real. After zero padding, the −N/2 modes have no +N/2 partner,
so g and h get imaginary parts. The packed product is then wrong at first order in those parts
(g_r·g_i − h_r·h_i), and this does not cancel when G ≠ H (the B1 and B2 terms).
On the synthetic test inputs, which have nonzero Nyquist modes, the `INFO` rows show the size of the error:
at N=32–64, equil and orthoCMB change by O(1) on non-Nyquist modes, and local by about 5%.
`KERNEL_ZERO_NYQUIST=1` removes it, and then M=3N/2 and M=2N agree to round-off.

**This does not affect real runs.** The Zeldovich_MPI dumps already have every Nyquist mode exactly zero
(checked on `1rank-run-cosm500/bin_files`). There, compat mode and `KERNEL_ZERO_NYQUIST=1` give identical output,
so M=3N/2 is safe for production. It only matters if Φ ever comes from a generator that keeps the Nyquist modes.
