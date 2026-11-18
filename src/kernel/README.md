# Kernel_MPI

`Kernel_MPI` computes the primordial non-Gaussianity (f_NL) kernels used to build non-Gaussian
initial conditions. It is step 2 of four:

1. **Dump pass** (`Zeldovich_MPI` built with `DUMP_GAUSSIAN_D`) writes the Gaussian potential as
   y-slices: `Phi_gaussian_y.{y}`.
2. **Kernel pass** (`Kernel_MPI`, this code) reads those slices and writes one set of kernel
   slices for each enabled f_NL shape: `kernel_<shape>_y.{y}`.
3. **Combine** (`PNG_Combine`) reads `Phi_gaussian` and the kernel slices and writes the white
   noise `(phi_G + f_NL * kernel) * s`.
4. **White-noise pass** (`Zeldovich_MPI` built with `LOAD_WHITE_NOISE_FROM_FILE`) builds the
   particles from that white noise.

[`1rank-run-cosm500/PPD256_fnl0_cluster.md`](../../1rank-run-cosm500/PPD256_fnl0_cluster.md)
walks through all four steps on a cluster. [`docs/png-pipeline.md`](../../docs/png-pipeline.md)
gives the full chain along with the normalization at each step.

## Build

`Kernel_MPI` is built by Meson along with `Zeldovich_MPI` (the target is in
[`src/meson.build`](../meson.build)). It links against `libzeldovich_ic` for the FFT, MPI exchange
and decomposition code.

The binary goes to `<build dir>/src/Kernel_MPI`, in whatever build directory you use. A full
`meson compile` builds it along with everything else, so
[`1rank-run-cosm500/build.sh`](../../1rank-run-cosm500/build.sh) already produces
`build/cosm500_1d_1rank/src/Kernel_MPI`. To build it on its own in a separate directory:

```bash
meson setup build/kernel                  # once; add --reconfigure if the dir exists
meson compile -C build/kernel Kernel_MPI  # -> build/kernel/src/Kernel_MPI
```

Compute precision follows `-Ddouble_precision`, which defaults to single. File precision is set
separately (see [Inputs and outputs](#inputs-and-outputs)), and values are converted on read and write.

## Run

```bash
mpiexec -n R build/kernel/src/Kernel_MPI [param.par] [--png-dir DIR] [--out-dir DIR] [--pad-2N]
```

The parameter file defaults to `.runtime/param.par`. It must be a plain `.par` file, because
ParseHeader does not evaluate `#include`, `@Var@` or expressions. Expand a `.par2` file first with
GenParam or `scripts/par2_to_par.py` (see `1rank-run-cosm500/run.sh`).

**Rank count:** `R` must be a multiple of `NumZRanks` and at most `N/2 + 1`. It does not have to
match the dump run or `PNG_Combine`, because files are keyed by y only.

### Parameter-file keys

| Key | Required | Meaning |
|---|---|---|
| `BoxSize` | yes | Box size; sets `kfund = 2*pi/BoxSize` |
| `NumZRanks` | yes | z-dimension of the rank grid (same Cartesian topology as `Zeldovich_MPI`) |
| `ZD_f_NL_local`, `ZD_f_NL_equil`, `ZD_f_NL_orthoCMB`, `ZD_f_NL_orthoLSS` | at least one nonzero | Every nonzero shape is computed, one after another |
| `ZD_PNG_DIR` | no (default `png_slices`) | Input directory, also the default output directory |
| `KernelPadM` | no | Padded FFT grid size `M` (see [Padding](#padding)) |
| `ZD_PNG_FileFloat64` | no (default `0`) | `1` writes the kernel files as float64, `0` as float32 |

Any other keys are ignored. Only whether an f_NL value is nonzero matters here: its value is
applied later, by `PNG_Combine --set`.

### Command-line options and environment variables

Command-line arguments override environment variables, which override the parameter file.

| Option / variable | Effect |
|---|---|
| `PARAM_FILE` | Parameter file, if none is given on the command line |
| `--png-dir DIR` / `PNG_DIR` | Override `ZD_PNG_DIR` |
| `--out-dir DIR` / `OUT_DIR` | Write outputs here instead of the input directory |
| `--pad-2N` / `KERNEL_PAD_2N=1` | Default `M = 2N`, as in `kernel_MPI.py` |
| `KERNEL_PAD_M=<M>` | Explicit `M`; overrides `KernelPadM` and the default |
| `KERNEL_ZERO_NYQUIST=1` | Zero the `kx`, `ky`, `kz = N/2` modes of phi before computing kernels |
| `KERNEL_ONLY_TERM=<A1\|C1\|...>` | Debug: run one term with coefficient 1 and write `kernel_<term>_y.{y}` |
| `KERNEL_FORCE_STAGED_A2A=1` | Test: force the point-to-point (MPI 3) Alltoallv path instead of `MPI_Alltoallv_c` |

## Inputs and outputs

| | Files | Notes |
|---|---|---|
| In | `<PNG_DIR>/Phi_gaussian_y.{y}`, y = 0 … N/2 | `N` is detected from the size of `Phi_gaussian_y.0` |
| Out | `<OUT_DIR>/kernel_<shape>_y.{y}`, y = 0 … N/2 | `<shape>` is `local`, `equil`, `orthoCMB` or `orthoLSS` |

Each file is one N×N complex slice of interleaved `(real, imag)` values, x fastest and then z, so
element `(x, z)` is at index `z*N + x`. Only the slices with y ≤ N/2 are stored; the other half
follows from Hermitian symmetry.

The element type is not stored in the file. Readers infer it from the file size: `N*N*8` bytes is
float32 and `N*N*16` bytes is float64 (`utils/png_slice_io.h`). Inputs may be either. Outputs are
float32 unless `ZD_PNG_FileFloat64 = 1`. Comparisons at the 1e-10 level need a double-precision
build *and* float64 files.

Rank 0 prints a run summary, the time for each term, a timing breakdown by stage, and peak memory.

## The method

Every term in `kernel_MPI.py` has the form

```
out += coef * J(k) * Q(k),   Q = FFT[ IFFT(G*phi) * IFFT(H*phi) ]
G = |k|^pG,   H = |k|^pH,   J = cJ * |k|^pJ
```

and each f_NL type is a weighted sum of these terms (from [`kernel_terms.c`](kernel_terms.c)):

| f_NL type | Terms |
|---|---|
| local | `A1` |
| equil | `-3 A1 - 6 C1 + 4 B1 + 2 B2` |
| orthoCMB | `-9 A1 - 24 C1 + 10 B1 + 8 B2` |
| orthoLSS | `A1, C1, B1, B2, D1, E2, F1, G1`, with coefficients derived for a scale-invariant spectrum |

The two inverse FFTs share a single complex transform: `z = (G + iH) * phi`, and then
`Q = Im(z^2) / 2`.

### Padding

The quadratic product is computed on a padded `M^3` grid to avoid aliasing. By default `M` is
the smallest valid size of at least `3N/2`, which is alias-free and smaller than the
`M = 2N` used by the Python code. `M` must satisfy `M >= N`, `M % 4 == 0`, and `M - N` even. If
`M < 3N/2` the program warns that the product will alias.

Padding and cropping are never done as separate steps. The original index `i` maps directly to
padded index `pad_idx(i)` (the same mapping as `fftshift -> np.pad -> ifftshift`).

### Parallel layout

There are two decompositions:

- **y-slices (original N grid):** the N/2+1 conjugate pairs `(y, N - y)` are split across ranks.
  Each rank reads, owns and writes its own `Phi_gaussian_y.{y}` / `kernel_<shape>_y.{y}`
  files. Mirror slices are built on the fly from Hermitian symmetry.
- **(x, z) pencils (padded M grid):** each rank owns a block of full-length y columns.

Each term runs in three stages ([`kernel_pass.h`](kernel_pass.h)):

1. **A, forward exchange:** fill `(G + iH) * phi` at the padded positions, run a 2D backward FFT
   on each y-slice, then Alltoallv into pencils.
2. **B, per z-slab:** 1D backward FFT in y, form `Q`, 1D forward FFT in y, then keep
   `y <= N/2`.
3. **C, return exchange:** Alltoallv back to the slice owners, run a 2D forward FFT, crop to N×N,
   then `out += coef * J * Q`.

To save memory, `q_half` reuses the persistent `recv_buffer`, and the shapes are computed one
after another into the same output buffers. Peak memory is therefore the same as for a
single-shape run.

## Conventions and normalization

These are the conventions the output depends on. `PNG_Combine` assumes all of them.

**Input field.** `Phi_gaussian_y.{y}` holds the Gaussian Bardeen potential written by the
`Zeldovich_MPI` dump pass, `φ(k) = (3/5) · sqrt(P_prim(k) / P_δ(k)) · δ(k)`, where `δ(k)` is the
Gaussian density mode (`generation/ZD_MPI_generation.c`). The generator's Fourier convention is
the unnormalized FFTW backward transform,

```
φ(x) = Σ_k φ(k) e^{+i k·x},    <|φ(k)|²> = P_Φ(k) / L³,
```

and the `1/L³` comes from `PowerSpectrum::Normalize` (`deps/zeldovich_core/src/power_spectrum.cpp`).
The DC mode, every Nyquist plane (`|kx|`, `|ky|` or `|kz| = N/2`) and modes past `ZD_k_cutoff` are
zero in real dumps. `KERNEL_ZERO_NYQUIST` therefore changes nothing for real dumps. It matters for
synthetic inputs, where the Nyquist planes are not zero.

**What a term computes.** With unnormalized FFTW transforms on the padded `M³` grid,

```
g(x) = Σ_k G(k) φ(k) e^{+ik·x},   h(x) = Σ_k H(k) φ(k) e^{+ik·x}      (backward FFTs)
Q(k) = (1/M³) Σ_x g(x) h(x) e^{-ik·x}                                 (forward FFT)
     = Σ_{k'} G(k') φ(k') H(k-k') φ(k-k')                             (alias-free if M ≥ 3N/2)
```

The `1/M³` cancels the padding, so `Q` does not depend on `M`. `Q` is in the same convention as
`φ`, so `g(x) h(x) = Σ_k Q(k) e^{+ik·x}` on the original grid. This is why `PNG_Combine` can form
`φ_G + f_NL · K` directly, with no extra factors. For `local` (`A1`, `G = H = J = 1`), `K` is the
transform of `φ(x)²`.

**Mean subtraction.** The `k = 0` mode of `Q` is set to zero (`stage_return`). For local this gives
`K = φ² − <φ²>_box`, which is the usual `Φ = φ + f_NL (φ² − <φ²>)`.

**Powers of k at k = 0.** `kpow(0, p)` is `1` for `p = 0` and `0` for every other `p`, including
negative ones (`kernel_terms.h`). This matches `kernel_MPI.py`. The DC input mode is zero anyway,
so this only stops `0 · ∞` from producing a NaN.

**Units.** `|k| = (2π / BoxSize) · sqrt(nx² + ny² + nz²)`, with `nx = kfreq(ix, N)` the signed
integer frequency. Only `BoxSize` enters, through `kfund`. Shapes that are homogeneous of degree 0
in `k` overall (`pG + pH + pJ = 0`, true of every shape in the table) do not depend on `BoxSize`.

**Shape coefficients.** The `equil`, `orthoCMB` and `orthoLSS` term coefficients are copied from
`kernel_MPI.py`. The `orthoLSS` coefficients assume a scale-invariant spectrum (`n_s = 1`)
(`kernel_terms.c`). 
We took these kernel shapes from https://arxiv.org/pdf/2206.01619 which uses the kernel generation approach from https://arxiv.org/pdf/1108.5512 

## Memory

The main per-rank allocations, in units of `sizeof(fftw_complex_t)` (8 bytes single, 16 double),
are:

| Buffer | Size | Lifetime |
|---|---|---|
| `phi`, `out` | `2 · pairs_per_rank · N²` | whole run |
| `recv_buffer` / `q_half` | `N · M² / R` (approximately) | whole run |
| `local_y_slices` + forward `send_buffer` | `4 · M²` | stage A |
| return `send_buffer` + `recv_buffer` + `q_plane` | about `3 · M²` | stage C |

The `M²` stage buffers do **not** shrink as ranks are added. At large `N` they set a floor on
memory per rank. For example, at `N = 12000` and `M = 18000`, single precision, stage A alone needs
about 10 GB per rank.

## Limitations

- `R ≤ N/2 + 1`: each rank must own at least one conjugate pair of y-slices.
- The `M²` per-rank buffers above set the largest `N` that fits on a node, however many ranks are used.
- The `orthoLSS` template is only exact for `n_s = 1`.
- Slice files have no header, so their type is inferred from their size. A truncated float64 file
  of exactly half the size would be read as float32.

## Source files

| File | Contents |
|---|---|
| [`kernel_main.c`](kernel_main.c) | Driver: arguments, parameters, topology, reading slices, the loop over shapes and terms, writing output, timing |
| [`kernel_params.cpp`](kernel_params.cpp) / `.h` | Reads the `.par` file with ParseHeader |
| [`kernel_grid.c`](kernel_grid.c) / `.h` | Choosing `M`, `pad_idx`/`kfreq`, pair ownership and pencil layout (`KernelLayout`) |
| [`kernel_pass.c`](kernel_pass.c) / `.h` | One kernel term: stages A–C, FFT plans, buffers (`KernelPass`) |
| [`kernel_terms.c`](kernel_terms.c) / `.h` | Shape table (`A1` … `G1`) and the f_NL term coefficients |
| [`kernel_io.c`](kernel_io.c) / `.h` | Slice file read/write, Hermitian fix-ups for the y=0 and y=N/2 planes, Nyquist zeroing |
