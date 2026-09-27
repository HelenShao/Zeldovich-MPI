# Loading Willem's (FLAMINGO / monofonIC-Panphasia) white noise for massive-neutrino runs

Notes for planning how the Zeldovich-MPI code reads external white noise and applies our own CLASS power spectrum. The goal is to generate ICs with the same phases as FLAMINGO for our massive-neutrino cosmology runs.

Related material:
- `CODE_REVIEW.md` Q17: the white-noise reader and the slab file format.
- `FLAMINGO_test.conf`: Willem's example config.
- `../monofonic`: a shallow clone of `wullm/monofonic@output_k_space_slabs`, commit `2de3faa`.
- `validate_external_noise.py`: the validation tests below.

## Changelog

| Date | Change |
|---|---|
| 2026-09-26 | Created. Moved `CODE_REVIEW.md` Q18 here. Answers on the half-cell phase, the `mode_weightings` normalization, zwindstroom, `KK1025` fixing, `fixed_power`, the Fix 4 guard, and the RNG draws. Three validation tests, plus `validate_external_noise.py` writing its numbers to a text file. |
| 2026-09-26 | New §3a: the N = 256 validation plan (no whitening step; recolour \(\sqrt{P_{\rm ours}}\,w\)), with one-line summaries of `fixed_power` and the Fix 4 guard. |
| 2026-09-26 | New §3b (why the slabs are white noise; \(N_0\) derivation, FLAMINGO \(N_0=2048\); `zeldovich_ps_inputted_ic` missing, implementation given; `fixed_power` is global) and §3c (build checklist). |
| 2026-09-26 | New §3d: `LOAD_EXTERNAL_NOISE` / `DUMP_D_SLABS` mode implemented on branch `PNG` (wrapper `zeldovich_ps_scaling`, slab reader with crop + transpose, in-loop mask + colouring, D dump), PNG code untouched. |
| 2026-09-26 | Colouring function renamed to `zeldovich_ps_scaling`; §3d notes on how it relates to `zeldovich_ps_inputted_ic` and on Gaussian-only runs. |
| 2026-09-26 | New §3d-bis: branch `neutrinos` off `meson` with only the external-noise mode (PNG WIP committed on `PNG`). |

---

## 1. What Willem's files are

- **Writer.** In Willem's branch, `external/panphasia_ho/pan_mpi_routines.c` (end of `PANPHASIA_compute_kspace_field_`, ~L520–610) writes Panphasia white noise straight from the random-field plugin, **before any transfer function or power spectrum**.
- **Files.** `output_k_space_slab_y.{0..N0/2}`: one file per \(k_y \ge 0\), each holding \(N_0^2\) complex values.
  - Values are raw `FFTW_REAL` pairs `(re, im)`: double by default (`CODE_PRECISION=DOUBLE`), float if monofonIC was built single-precision.
  - Order is `[x][z]` with z fastest. There is no header and no labels.
- **Hermitian symmetry.** The \(k_y=0\) plane is already Hermitian (z > N/2 filled with conjugates). Nyquist modes are zero, so the last slab is all zeros.
- **Grid size.** \(N_0 = f_{\rm dim}\cdot(\text{descriptor base size}\ll\text{level})\), with `fdim = 2` in `src/plugins/random_panphasia.cc` L136.
  - `N0` can be larger than `GridRes`. monofonIC then Fourier-interpolates down to `GridRes` (`FourierInterpolateCopyTo`).
  - The monofonIC log prints `N0`.
- **L1000N1800.** FLAMINGO naming: a 1000 Mpc box (673.2 Mpc/h, \(h=0.6732\)) with \(1800^3\) DM particles, PLANCK \(\Sigma m_\nu = 0.06\) eV.
  - The "L1000N1800 slab files" are what `FLAMINGO_test.conf` produces with this branch.
  - They contain exactly the white noise we need.

## 2. Q&A

### 2.1 The half-cell phase factor: should we apply the same convention?

**What it is** (`pan_mpi_routines.c` L384–387). Every mode is multiplied by
\[
e^{-i\pi f_{\rm dim}(k_x+k_y+k_z)/N_0}=e^{-i\mathbf k\cdot\boldsymbol\Delta},\qquad \boldsymbol\Delta=\tfrac{f_{\rm dim}}{2}\tfrac{L}{N_0}(1,1,1).
\]
This translates the real-space field by half a Panphasia coefficient cell, which is one cell of the \(N_0\) FFT grid since `fdim = 2`. It places the field at the centres of Panphasia's cells rather than their corners. It changes neither \(|w|\) nor the Hermitian symmetry, so P(k) and all other statistics are unchanged.

**Should we multiply by \(e^{+i\mathbf k\cdot\boldsymbol\Delta}\)? No: use the loaded \(w\) exactly as it is.** The phase factor is part of the Panphasia field that monofonIC itself uses. monofonIC then does two things:
- Uses \(w\) as-is on its grid (after cropping to `GridRes`).
- Puts SC-lattice particle \((i,j,k)\) at \((i,j,k)/N\) with zero shift (`include/particle_generator.hh` L41, `grid_fft.hh` `get_unit_r`). So particle \((i,j,k)\) gets the displacement sampled at grid point \((i,j,k)\).

Our code does the same: particle \((i,j,k)\) gets grid point \((i,j,k)\)'s displacement. Taking \(w\) unchanged therefore reproduces monofonIC's particle-by-particle correspondence, up to our LPT order and P(k).

What remains can differ by a *global* coordinate origin, e.g. Abacus's box-centred \([-L/2, L/2)\) vs monofonIC's \([0, L)\). That is a rigid translation of the whole periodic box and is handled when comparing positions. Test 3 checks for a residual shift; the expected answer is 0 cells.

### 2.2 The `/ sqrt(ptr_mode_weightings[index1])` placement

The line is
```c
phase_shift_and_scale = sqrt(fdim^3) * cexp(fdim*(-I)*pi*(kx+ky+kz) / sqrt(ptr_mode_weightings[index1]) / nfft_dim) / pow(nfft_dim, 1.5);
```
As written, the division by \(\sqrt{\text{weights}}\) is **inside** the exponent. That has two effects:
1. The amplitude is *not* renormalized by the mode weights.
2. The phase picks up a small mode-dependent extra term, proportional to \(1/\sqrt{\text{weights}}-1\).

The weights are the summed \(|\text{spherical-Bessel weight}|^2\) of the Panphasia basis. They are about 1 for well-resolved modes and drop towards Nyquist, so both effects are small and grow near Nyquist. If the division had been meant as an amplitude normalization, it would sit outside `cexp()`.

**The identical line is in upstream monofonIC** (`bitbucket.org/ohahn/monofonic`, `external/panphasia_ho/pan_mpi_routines.c` L383–386), so Willem didn't introduce it.

Consequences:
- **Matching FLAMINGO:** use the files as-is, because this is what FLAMINGO's ICs used.
- **Our recolouring:** the normalization matters. Test 1 measures \(\langle|w|^2\rangle\) overall and per \(k\)-shell, which would show any drop near Nyquist.
- **Ask Willem** whether it's intentional.

### 2.3 Do we need zwindstroom to get the Panphasia white-noise files?

**No.** `transfer = zwindstroom` is monofonIC's transfer-function plugin. It uses:
- **CLASS**, for the transfer functions;
- **zwindstroom** (Willem's code), for the scale-dependent growth of a 3-fluid model (baryons, CDM, neutrinos) in N-body gauge.

In this branch they are downloaded by CMake `FetchContent` at configure time (`external/class.cmake`, `external/zwindstroom.cmake`), not git submodules, which needs internet on the build node. The slab files are written inside the random-number plugin, before any transfer function runs. To generate the files, use `transfer = eisenstein` and build with `-DENABLE_CLASS=OFF -DENABLE_ZWINDSTROOM=OFF`. The white noise is identical, because it depends only on the Panphasia descriptor and the grid.

Relevant later: what zwindstroom computes (scale-dependent growth \(D(k,z)\) and growth rate \(f(k,z)\) with neutrinos, N-body gauge) is exactly what our own P(k) and velocities must handle (§4).

### 2.4 Why is \(|w| = 1\) for modes with \(k^2 \le 1025\)?

The descriptor contains `KK1025`, which sets Panphasia's `descriptor_kk_limit`. For every mode with integer \(k_x^2+k_y^2+k_z^2 \le 1025\) (i.e. \(|\mathbf k| \lesssim 32\) fundamental modes, excluding \(k=0\) and Nyquist), the code divides by \(|w|\) (`pan_mpi_routines.c` ~L401–428). So those modes keep their random phase but have unit amplitude.

This is **partial amplitude fixing** (Angulo & Pontzen 2016). On the largest scales there are few modes per \(k\)-bin, so their random amplitudes cause large cosmic variance in P(k) and in anything sensitive to large scales (e.g. BAO, large-scale bias). Fixing the amplitudes removes that scatter. Limiting it to \(k^2 \le 1025\) keeps the field Gaussian on smaller scales, where fixing could bias non-linear statistics and higher-order statistics. This is FLAMINGO's choice for its ICs.

### 2.5 Why not \(|D|=\sqrt{P}\) for every mode (i.e. `fixed_power = 1`)?

1. **To match FLAMINGO.** FLAMINGO fixed only \(k^2 \le 1025\). The file already encodes that: those modes have \(|w|=1\), and the rest are Gaussian. With `fixed_power = 1`, our code would force \(|D| = \sqrt P\) on *all* modes, and our ICs would no longer match FLAMINGO's small-scale amplitudes.
2. **Physics.** A fully fixed field isn't Gaussian. Mean P(k) is fine, but covariances, the one-point PDF and higher-order statistics are affected, and non-linear mode coupling on small scales changes the evolved field slightly. Fixing is a variance-reduction trick best kept to scales where it's harmless.

So use `fixed_power = 0` and take \(|w|\) from the file. Also check that `zeldovich_ps_inputted_ic` doesn't apply `fixed_power` itself.

### 2.6 Why is the Fix 4 mask guard needed if Nyquist modes are already zero in the file?

Several of our masked modes are *not* zero in the file:
1. **Cropping from \(N_0\) to \(N\).** The file's zeros sit at *its* Nyquist, \(N_0/2\). Our Nyquist is \(N/2\) (900 for N = 1800). If \(N_0 > N\), the file has non-zero values at \(|k_i| = N/2\), and those must be zeroed.
2. **Spherical cutoff.** With `CornerModes = 0`, our code zeroes \(k^2 \ge (N/2)^2/k_{\rm cut}^2\) (the corners of the cube). The file doesn't. To *match FLAMINGO*, which keeps the corners, use `CornerModes = 1`, `k_cutoff = 1`.
3. **The DC mode**, which should be exactly 0 in any case.
4. **PNG runs.** The \(f_{\rm NL}\) kernel is non-zero at Nyquist and cutoff modes no matter what the file contains.
5. **Robustness.** The generator shouldn't depend on any particular input file being clean.

The guard costs nothing and makes masked modes identical to the Gaussian path.

### 2.7 How are the RNG draws "identical" if we import external white noise?

They aren't used at all. In the load path the code still calls `get_cgauss` and advances `nskip` for every mode, but the drawn value is thrown away and replaced by the file value. The draws are kept only because the loop was left unchanged.

"Identical random stream" mattered only for Grace's `LOAD_D_FROM_FILE` round trip, where \(\phi_G\) had to be regenerated bit-for-bit to add \(f_{\rm NL}K\). For external noise it's irrelevant. A cleaner external-noise path would skip `get_cgauss` entirely (faster) and drop the dependence on `ZD_Seed` and `MAX_PPD`.

---

## 3. Plan: reading external noise and applying our own CLASS P(k)

1. **Reader** (replaces Grace's labelled-record reader):
   - read `output_k_space_slab_y.{y}` for \(y \le N/2\), with `fread` of \(2N_0^2\) `FFTW_REAL`;
   - map `(ix, iz)` → `[iz*N + ix]`;
   - if \(N_0 > N\), crop by signed index (target \(k \in (-N/2, N/2)\); source index \(k\ge0\ ?\ k : k+N_0\));
   - call `MPI_Abort` on a missing or short file;
   - files are keyed by y, so any rank count works.
2. **Recolour** with \(D=\sqrt{P_{\rm code}(k)}\,w\). \(P_{\rm code}\) is `zeldovich_ps_power(ps_handle, k)`, the same quantity `get_cgauss` uses. `ps_handle` already carries our CLASS P(k) file, cosmology, normalization and growth to \(z_{\rm start}\).
3. **Parameters to match FLAMINGO:**
   - `N = 1800` (or the target), `BoxSize = 673.2` Mpc/h;
   - `CornerModes = 1`, `k_cutoff = 1`, `fixed_power = 0`;
   - the same \(z_{\rm start}\) (31) if comparing directly.
4. **Unchanged:**
   - the Fix 4 mask guard;
   - the self-conjugate mirroring (the file is already Hermitian and uses the same "keep z < N/2" half, so this rewrites identical values);
   - displacements, PLT, FFTs, MPI and output.
5. **Optional:** skip the `get_cgauss` draws in this mode (§2.7). Add an input switch, e.g. `ZD_ExternalNoiseMode = white | delta`, in case a file already holds \(\delta\).
6. **PNG on external noise:** colour with \(P_{\rm prim}\), apply the \(f_{\rm NL}\) kernel, whiten, then recolour with \(P_\delta\). The normalization from §2.2 matters here.

### 3a. Proposed small validation run (N0 = PPD = N = 256)

1. **Run Willem's branch** with `GridRes = 256` and the FLAMINGO descriptor.
   - Confirm in the log that `N0 = 256` (\(N_0 = f_{\rm dim}\cdot2^{\rm level}\), `fdim = 2`).
   - Transpose `[x][z]` → `[z][x]` (x fastest), preferably inside the C reader rather than by rewriting files.
   - Panphasia phases don't depend on resolution, so this run shares FLAMINGO-1800's large-scale modes.
2. **No whitening step.** The slab files are already unit white noise, written inside Panphasia before any transfer function. Dividing by a FLAMINGO \(\sqrt P\) would corrupt them.
   - Grace's `LOAD_WHITE_NOISE_FROM_FILE` also doesn't whiten: it reads \(w\) and applies \(\sqrt{P_\delta}\).
   - Whitening (\(\times\tfrac32/\sqrt{P_{\rm prim}}\)) exists only in the PNG `LOAD_D_FROM_FILE` path.
   - It's only needed if a file holds \(\delta\) (the `ZD_ExternalNoiseMode = delta` case).
3. **Recolour** with \(D=\sqrt{P_{\rm ours}(k)}\,w\), using our CLASS P(k). This matches `get_cgauss` (`zeldovich_wrapper.cpp` L303–317): \(D=\sqrt P\times\)(unit complex Gaussian).
4. **Plug `D` into the unchanged pipeline**, with `CornerModes = 1`, `k_cutoff = 1`, `fixed_power = 0` and the Fix 4 guard. Run tests 1–3, saving the numbers to a text file with `--out`, then the Abacus simulation.

**`fixed_power` in one line.** `0` gives \(D=\sqrt P\,w\), which keeps FLAMINGO's mixture (\(|w|=1\) for \(k^2\le1025\), random elsewhere). `1` forces \(|D|=\sqrt P\) on every mode, erasing the random small-scale amplitudes. Ensure `zeldovich_ps_inputted_ic` ignores `fixed_power`.

**Fix 4 guard in one line.** Modes the Gaussian path zeroes (DC, \(|k_i|=N/2\), and the spherical cutoff if `CornerModes = 0`) stay zero: `if (!is_masked) D = sqrt(P)*w_loaded;`. At N0 = N = 256 with `CornerModes = 1` it only matters for DC. It's essential when cropping (N0 > N), when `CornerModes = 0`, and for PNG.

### 3b. Clarifications (2026-09-26)

- **Why the slabs are white noise.**
  - The writer is inside `PANPHASIA_compute_kspace_field_`, called from `RNG_panphasia::Fill_Grid`.
  - In `src/ic_generator.cc`, `Fill_Grid(wnoise)` (L279) runs before monofonIC's normalization/fixing (`wn / volfac`, L361–367) and before any transfer function.
  - Unit variance is strongly indicated (the `KK` block sets \(|w|=1\); the \(\sqrt{f_{\rm dim}^3}/N_0^{1.5}\) scale). It is confirmed by test 1, including per-shell \(\langle|w|^2\rangle\) because of the `mode_weightings` quirk.
- **\(N_0\).** The descriptor `S1` means the box is 1 cell wide at level `L18`, so `descriptor_base_size = 1` (`high_order_panphasia_routines.c` L1289, L1334).
  - At `rel_level` levels deeper the coefficient grid is \(2^{\rm rel\_level}\) per side, and the FFT grid is `fdim = 2` times finer, so \(N_0 = 2\cdot2^{\rm rel\_level}\).
  - It is the smallest value \(\ge\) `GridRes` (`random_panphasia.cc` L146).
  - `GridRes = 256` gives \(N_0 = 256\) (no crop). `GridRes = 1800` gives \(N_0 = 2048\): FLAMINGO slabs are 1025 files of \(2048^2\) complex doubles (~69 GB), cropped to \(|k_i|<900\).
- **Transpose:** in the C reader, `(ix, iz)` → `[iz*N + ix]`.
- **`zeldovich_ps_inputted_ic` does not exist** in the repo or Grace's folder; it is only called. Implement it in `utils/zeldovich_wrapper.{h,cpp}`:
  ```cpp
  void zeldovich_ps_inputted_ic(PowerSpectrumHandle ps, double wavenumber,
                                double w_re, double w_im, double* re, double* im) {
      if (!ps || !re || !im) return;
      PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
      const double A = sqrt(p->power(wavenumber));   // same P as zeldovich_ps_cgauss_from_buffer
      *re = A * w_re;
      *im = A * w_im;
  }
  ```
  It has no `fixed_power` branch: \(|w|\) comes from the file.
- **`fixed_power` is global** (`p->fixed_power`, `zeldovich_wrapper.cpp` L310): all modes fixed or none. FLAMINGO's per-mode selection comes from Panphasia's `KK1025` inside the white noise (monofonIC `DoFixing = no`), so with external noise it is inherited from the file.

### 3c. Build checklist before the N = 256 run

1. Implement `zeldovich_ps_inputted_ic` (§3b) plus its header declaration.
2. Reader for `output_k_space_slab_y.{y}`:
   - configurable directory and precision;
   - transpose to `[z][x]`;
   - crop for \(N_0>N\);
   - `MPI_Abort` on a missing or short file.
3. `config.h`: add `LOAD_WHITE_NOISE_FROM_FILE`. Resolve the `ps_handle_primordial_dimensional` signature: wire it through the driver, or put it behind a PNG-only `#if` (the external-noise path doesn't need it).
4. Fix 4 mask guard in the load block.
5. Dump the recoloured `D` (y-keyed, `[x][z]`, complex double) for test 2.
6. Run parameters:
   - `N = 256`, `BoxSize = 673.2`;
   - `CornerModes = 1`, `k_cutoff = 1`, `fixed_power = 0`;
   - our \(P_{cb}(k,z_{\rm start})\) file;
   - check the P(k) rescaling options (\(\sigma_8\), growth), since `p->power(k)` is used as-is.
7. monofonIC: one run with `GridRes = 256`, the FLAMINGO descriptor and `[testing] test = potentials_and_densities`. It writes the slabs (during `Fill_Grid`) and the first-order `delta` HDF5 for test 3. Needs Panphasia licence registration and an MPI + FFTW3-MPI build.
8. Run `validate_external_noise.py --out run256/validation_results.txt`, then the Abacus run.

### 3d-bis. Branch `neutrinos` (2026-09-26) — use this one for the neutrino runs

- `PNG` work committed locally as `c9263ae` (WIP: Grace's blocks + external-noise mode). Not pushed.
- `neutrinos` branched from `meson` (`c46448a`); commit `60834ac` adds only the external-noise mode, no PNG code. It builds and links independently of Grace's missing pieces.
- Difference from the `PNG` version: the colouring sits in the existing unmasked `else` branch (the `meson` mask already zeroes DC, Nyquist and cutoff modes), so no duplicate mask check and no wasted `get_cgauss` draw; `get_cgauss` is marked `[[maybe_unused]]`.
- Verified locally: `ZD_MPI_generation.cxx` compiles with no warnings with the flags off and with `-DLOAD_EXTERNAL_NOISE=1 -DDUMP_D_SLABS=1` (`-Dextra_cpp_args=-DLOAD_EXTERNAL_NOISE=1,-DDUMP_D_SLABS=1`). A full local link was not possible: the ParseHeader subproject's bison step doesn't produce `phParser.tab.cc` on this Mac (toolchain issue, unrelated; build on Aurora).
- Reconciling later: merge `neutrinos` into `PNG`; the conflict will be in the in-loop block, where the `PNG` version needs its own mask check because Grace's load blocks overwrite D after the mask.

### 3d. Implemented on branch `PNG` (2026-09-26)

Items 1, 2, 4 and 5 of §3c are implemented as a separate `LOAD_EXTERNAL_NOISE` mode. Grace's PNG blocks and `zeldovich_ps_inputted_ic` calls are untouched; the new mode is exclusive with them (`#error`).

| Piece | Location |
|---|---|
| `zeldovich_ps_scaling(ps, k, w_re, w_im, &re, &im)`: \(D=\sqrt{P(k)}\,w\), ignores `fixed_power` | `src/utils/zeldovich_wrapper.{h,cpp}` (after `zeldovich_ps_cgauss_from_buffer`) |
| Flags `LOAD_EXTERNAL_NOISE`, `EXTERNAL_NOISE_DIR`, `EXTERNAL_NOISE_PREFIX`, `EXTERNAL_NOISE_FLOAT`, `DUMP_D_SLABS`, `DUMP_D_DIR` | `src/config.h`, new "EXTERNAL WHITE NOISE" block |
| `load_external_noise_slab(global_y, N, w_re, w_im)`: infers \(N_0\) from the file size (must be even, \(\ge N\)), crops by signed wavenumber, transposes `[x][z]` to `[z*N + x]`, `MPI_Abort` on errors | `ZD_MPI_generation.c`, static function before `generate_zd_mpi_slice_pair_local` |
| Slab read before the z-loop | after Grace's `LOAD_WHITE_NOISE_FROM_FILE` reader |
| In-loop colouring with its own mask check (DC, Nyquist, cutoff sphere if `CornerModes = 0`) | after Grace's in-loop load blocks, before STEP 3 (F, G, H) |
| D dump `DUMP_D_DIR/D_slab_y.<y>`, `[x][z]`, complex double (the layout `validate_external_noise.py` expects) and frees | cleanup section, before `#undef PRIM_SLICE` |

Usage: compile with `-DLOAD_EXTERNAL_NOISE=1 -DDUMP_D_SLABS=1` (or edit `config.h`), put the slabs in `external_noise/` (or set `EXTERNAL_NOISE_DIR`), and run with `CornerModes = 1`, `k_cutoff = 1`. Each rank reads only the slabs of the y-values it generates. On the y = 0 and y = N/2 planes the dumped D is the raw coloured file values; the pipeline then overwrites half of those planes by Hermitian mirroring, which is consistent because the monofonIC y = 0 plane is already Hermitian.

Syntax checks (`mpicxx -std=c++17 -fsyntax-only`) pass with the flags at 0, with `LOAD_EXTERNAL_NOISE=1 DUMP_D_SLABS=1`, and with `PARALLELIZE_Z_LOOP=1`. Combining with `LOAD_D_FROM_FILE` triggers the `#error`.

`zeldovich_ps_scaling` vs Grace's `zeldovich_ps_inputted_ic`: same argument list and presumably the same intended math (\(\sqrt{P(k)}\times\) input), but hers has no definition anywhere, and her call sites apply no mask. Hers could be defined as a one-line forward to `zeldovich_ps_scaling`.

Gaussian-only runs (all PNG and external-noise flags 0): D is the usual `get_cgauss` draw and is unchanged, but the always-on `phi_G` block (two spline evaluations per mode, result unused) costs a little time, and the link fails until the signature below is fixed (e.g. pass `NULL` for the primordial handle).

Still blocking a full build (PNG side, deliberately left alone): the driver call does not pass `ps_handle_primordial_dimensional`, the header lacks it, `zeldovich_ps_inputted_ic` is undefined (only reached when `LOAD_WHITE_NOISE_FROM_FILE`/`LOAD_D_FROM_FILE` are on), and `config.h` lacks Grace's flags. With her flags off, the undefined function is compiled out, but the signature mismatch between definition, header and driver call will fail at link time.

## 4. Neutrino-specific items to settle in the plan

- **Which spectrum.** For CDM+baryon particles, use \(P_{cb}(k, z_{\rm start})\), not total matter \(P_m\). FLAMINGO's config has `WithNeutrinos = yes` with neutrino particles off. Decide between:
  - direct CLASS output at \(z_{\rm start}\) (synchronous gauge, or converted to N-body gauge);
  - back-scaling from \(z=0\) with scale-dependent growth, as monofonIC and zwindstroom do.
- **Velocities.** With massive neutrinos, the growth rate \(f(k)\) depends on scale. Zeldovich velocities \(\mathbf v \propto f\,\boldsymbol\psi\) with a single \(f\) are off at the per-cent level on large scales. Check what `Zeldovich-MPI` and PLT assume, and whether a \(k\)-dependent \(f(k)\) can be applied per mode.
- **LPT order.** FLAMINGO used 3LPT at \(z = 31\); we use 1LPT (+PLT). The phases match but the particles won't match exactly. Decide on the comparison metric (P(k) ratios, cross-correlation, halo matching).
- **Neutrino particles / \(\delta f\).** If our runs include neutrino particles, they need their own ICs, e.g. from the same white noise with the neutrino transfer function (FastDF in monofonIC). Out of scope for the CDM ICs, but the phases should be shared.
- **Units.** Mpc/h vs Mpc, and the \(h\) used in the CLASS file vs `BoxSize`.

## 5. Validation tests (run on a small grid first)

Use the same descriptor with a small `GridRes` (e.g. 128) so monofonIC's `N0` is small, and our code at the same \(N\). Script: `validate_external_noise.py`. It prints every number and **writes them to a text file** (`--out`, default `validation_results.txt`).

| Test | What it checks | Pass criteria |
|---|---|---|
| **1. White-noise statistics** (slab files only) | normalization \(\langle|w|^2\rangle\) overall and per \(k\)-shell; \(\mathrm{var(Re)}\), \(\mathrm{var(Im)}\); fixed modes; Nyquist; Hermitian \(y=0\) plane | \(\langle|w|^2\rangle\approx1\) (no drop near Nyquist, cf. §2.2); var ≈ 0.5 each; \(||w|-1| \approx 0\) for \(k^2\le1025\); Nyquist = 0; Hermitian violation ≈ 0 |
| **2. \(P_{\rm ours}(k)/P_{\rm input}(k)\)** (our linear \(D\) and the CLASS P(k) at \(z_{\rm start}\)) | recolouring normalization and units | shell-averaged ratio ≈ 1 (cosmic-variance scatter only); **exactly 1 for fixed modes** |
| **3. Cross-correlation with monofonIC's \(\delta\)** | same phases, axes, FFT sign, half-cell convention | \(r(k)\approx1\) at all \(k\); best-fit translation ≈ (0, 0, 0) cells; residual phase ≈ 0 |

**Inputs:**
- **Test 1:** `--slab-dir`, `--N0`, `--precision`, `--kk 1025`.
- **Test 2:** our \(D\) as slabs in the same layout (`D_slab_y.{y}`, `[x][z]`, complex double). This needs a small dump option in `ZD_MPI_generation.c` that writes `D` right after the load/recolour step, keyed by y. Also `--pk` (CLASS P(k) at \(z_{\rm start}\): \(k\) [h/Mpc], \(P\) [(Mpc/h)\(^3\)]) and `--box`.
  - Normalization: by default `P = V |D|^2`, i.e. \(\delta(\mathbf x)=\sum_k D\,e^{i\mathbf k\cdot\mathbf x}\), the unnormalized backward FFT. If our convention differs, pass `--ours-norm`; test 2 itself tells you whether the convention is right.
- **Test 3:** monofonIC run with `[testing] test = potentials_and_densities`, which writes first-order `delta` to `[output] fname_hdf5` (`src/testing.cc` L56–123). Pass `--mono-hdf5`. Only \(r(k)\) and the phase are compared, so normalization doesn't matter.
- **Diagnosing test 3 failures:**
  - an \(r(k)\) that falls like \(\cos(\mathbf k\cdot\boldsymbol\Delta)\), or a non-zero best-fit translation, points to an origin or half-cell convention mismatch;
  - \(r \approx 0\) points to swapped axes or conjugation (FFT-sign mismatch).

**Example** (keep the output file with the run):
```bash
python validate_external_noise.py \
  --slab-dir run128/slabs --N0 128 --kk 1025 \
  --ours-slabs run128/ours --N 128 --box 673.2 --pk run128/class_pk_cb_z31.dat \
  --mono-hdf5 run128/output.hdf5 \
  --out run128/validation_results.txt
```

The script was checked on synthetic data (N = 32): tests 1–3 pass, and an injected 1-cell shift is recovered as (+1, 0, 0) cells.

## 6. Questions for Willem

1. Can you send the L1000N1800 slab files (e.g. via Globus)? What are `N0` and the precision?
2. Is the `/sqrt(ptr_mode_weightings)` inside `cexp()` intentional (§2.2)? What \(\langle|w|^2\rangle\) should we expect near Nyquist?
3. For FLAMINGO, did monofonIC zero the \(|k_i| = N/2\) modes after cropping from `N0` to 1800?
4. Which spectrum and gauge did the FLAMINGO DMO ICs use for CDM (\(P_{cb}\), N-body gauge, back-scaled from \(z=0\) with zwindstroom)? Can monofonIC dump that \(P(k, z_{\rm start})\) so we can use the same one?
