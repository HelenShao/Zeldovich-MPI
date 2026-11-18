# Primordial non-Gaussianity pipeline

Non-Gaussian initial conditions are made in four steps by three programs. Each step reads and
writes y-slice files in `ZD_PNG_DIR`, so steps can use different rank counts.

```
 Zeldovich_MPI        Kernel_MPI                  PNG_Combine                         Zeldovich_MPI
 (dump build)         (src/kernel)                (src/png_combine)                   (white-noise build)
 Phi_gaussian_y.{y} ──▶ kernel_<shape>_y.{y} ──▶ phi_<shape>_fNL<f>_divided_y.{y} ──▶ particles
 φ_G                     K                        w
```

| Step | Program and build | Reads | Writes |
|---|---|---|---|
| 1. Dump | `Zeldovich_MPI`, `DUMP_GAUSSIAN_D=1` (the default) | P(k), P_prim(k) | `Phi_gaussian_y.{y}` (and Gaussian particles) |
| 2. Kernels | `Kernel_MPI` ([README](../src/kernel/README.md)) | `Phi_gaussian_y.{y}` | `kernel_<shape>_y.{y}` for every nonzero `ZD_f_NL_<shape>` |
| 3. Combine | `PNG_Combine` | φ_G, K, P_prim | `phi_<shape>_fNL<f>_divided_y.{y}` |
| 4. Particles | `Zeldovich_MPI`, `LOAD_WHITE_NOISE_FROM_FILE=1` | the step 3 file | particle ICs |

The two `Zeldovich_MPI` modes are compile-time macros, and only one can be on in a build. Give
each build its own directory. `Kernel_MPI` and `PNG_Combine` are built in either one:

```bash
meson setup build/dump                                                    # DUMP_GAUSSIAN_D=1 by default
meson setup build/noise -Dextra_cpp_args="-DDUMP_GAUSSIAN_D=0,-DLOAD_WHITE_NOISE_FROM_FILE=1"
```

[`1rank-run-cosm500/PPD256_fnl0_cluster.md`](../1rank-run-cosm500/PPD256_fnl0_cluster.md) runs
the whole chain on a cluster and checks the `f_NL = 0` round trip.

## The math, step by step

All fields follow the [Zeldovich_MPI conventions]: unnormalized FFTW
backward transforms (`f(x) = Σ_k f(k) e^{+ik·x}`), and every power spectrum scaled by `1/L³`.
`n(k)` is the unit complex Gaussian the generator draws, `<|n(k)|²> = 1`.

**1. Dump.** The generator draws `δ(k) = sqrt(P_δ(k)) n(k)` as usual and writes

```
φ_G(k) = (3/5) · sqrt(P_prim(k) / P_δ(k)) · δ(k)  =  (3/5) · sqrt(P_prim(k)) · n(k)
```

This is the Gaussian Bardeen potential `Φ = (3/5) ζ`, with `P_prim = P_ζ` read from
`ZD_Pk_primordial_filename`. Masked modes (DC, Nyquist planes, beyond `ZD_k_cutoff`) are zero.

**2. Kernels.** For each enabled shape, `Kernel_MPI` computes a quadratic kernel `K[φ_G]` in the
same Fourier convention as `φ_G`, with its mean removed. For `local`, `K = φ_G² − <φ_G²>`. See
[the kernel conventions](../src/kernel/README.md#conventions-and-normalization).

**3. Combine.** The non-Gaussian potential is divided back down to white noise:

```
w(k) = s(k) · [ φ_G(k) + f_NL · K(k) ],    s(k) = (5/3) / sqrt(P_prim(k)),   s = 0 on masked modes
```

At `f_NL = 0` this gives back `w = n`, up to the precision of the slice files. The mask matters because `K` has power at every
mode, including the Nyquist planes and the corners. `s = 0` removes those modes so that the final
field has the same mode set as a Gaussian run.

**4. Particles.** The white-noise build reads `w` and sets `δ(k) = sqrt(P_δ(k)) · w(k)`. After that
it runs the usual Zel'dovich/PLT path. It still draws (and discards) the Gaussian random numbers,
so the RNG streams stay the same.

### What this implies

- **`P_prim` sets the strength of the non-Gaussianity. `P_δ` sets the matter amplitude.** `n`
  cancels out of the ratio, so `ZD_Pk_sigma` renormalization of `P_δ` has no effect on `φ_G`.
  `P_prim` must carry the physical `A_s` amplitude, in units of `length³` that match `BoxSize`
  (it is converted by `ZD_Pk_scale` like the matter P(k), and normalized only by `1/L³`, never by
  σ(R)).
- `P_δ` and `P_prim` should describe the same cosmology. Otherwise the implied transfer function
  `sqrt(P_δ / P_prim)` is wrong.
- `ZD_f_NL` is the amplitude. `ZD_f_NL_<shape>` only switches a shape on or off; its value does
  not matter.
- Sign convention: the code uses `Φ = +(3/5) ζ` and `Φ = φ + f_NL (φ² − <φ²>)`.


## Files

Every slice file is `N × N` complex values, interleaved `(re, im)`, x fastest then z (element
`(x, z)` at `z*N + x`), for `y = 0 … N/2`. The other half follows from Hermitian symmetry. There is
no header; float32 (`N*N*8` bytes) and float64 (`N*N*16` bytes) are told apart by size. Writers use
float32 unless `ZD_PNG_FileFloat64 = 1` (`src/utils/png_slice_io.h`).

`<f>` in the step 3 name is the f_NL given to `PNG_Combine --set` (and the `ZD_f_NL` that step 4
reads), written as an integer when it is whole and with `%g`
otherwise, e.g. `phi_local_fNL10_divided_y.3` or `phi_equil_fNL-2.5_divided_y.3`.

Disk use: each file set is about `N³ · 4` bytes in float32 (about 7 TB at `N = 12000`). Step 1
writes one set, step 2 writes one per shape, and step 3 writes one per (shape, f_NL) pair.

## Parameters

| Key | Used by | Meaning |
|---|---|---|
| `ZD_PNG_DIR` | all | Directory for the slice files (default `png_slices`) |
| `ZD_PNG_FileFloat64` | all writers | `1` writes float64 slice files |
| `ZD_Pk_primordial_filename` | 1, 3 | `P_ζ(k)` table (default `class_pk_primordial_dimensional.dat`). If it cannot be read, `P_prim = 0` and the run continues with a warning |
| `ZD_f_NL_local`, `_equil`, `_orthoCMB`, `_orthoLSS` | 2, 4 | Shape switches (nonzero = on). Step 4 needs exactly one |
| `ZD_f_NL` | 4 | Amplitude; it forms the file name that step 4 reads, so it must be one of the step 3 values |
| `ZD_white_noise_name` | 4 | Read `{name}_y.{y}` instead; only allowed when every shape is 0 |

`PNG_Combine` ignores `ZD_f_NL` and the shape switches. It takes its shapes and amplitudes on the
command line instead, and can produce several in one run:

```bash
mpiexec -n R build/src/PNG_Combine run.par --set local:10,-10 --set equil:150 [--out-dir DIR] [--overwrite]
```

It does not overwrite existing outputs unless `--overwrite` is given.
