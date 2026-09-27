#!/usr/bin/env python3
"""Validation tests for loading external (monofonIC/Panphasia) white noise.

Test 1  white-noise statistics of the slab files (normalization, fixed modes, Nyquist, Hermitian y=0 plane)
Test 2  P_ours(k) / P_input(k) from our linear density field (needs --ours-slabs and --pk)
Test 3  cross-correlation r(k) and phase offset between our delta and monofonIC's delta (needs --mono-hdf5)

All numbers are printed and written to --out (default: validation_results.txt).

Slab format (Willem's output_k_space_slabs branch and, for our dump, the same layout):
  one file per k_y = 0..N/2, raw complex values (re, im) in [x][z] order with z fastest, no header.

Example:
  python validate_external_noise.py --slab-dir slabs/ --N0 128 --kk 1025 \
      --ours-slabs ours/ --N 128 --box 673.2 --pk class_pk_zstart.dat \
      --mono-hdf5 output.hdf5 --out validation_results.txt
"""
import argparse
import os
import sys

import numpy as np


class Tee:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *args):
        line = " ".join(str(a) for a in args)
        print(line)
        self.f.write(line + "\n")
        self.f.flush()


def read_slab(path, N, dtype):
    raw = np.fromfile(path, dtype=dtype)
    if raw.size != 2 * N * N:
        raise ValueError(f"{path}: expected {2*N*N} values, found {raw.size}")
    return (raw[0::2] + 1j * raw[1::2]).reshape(N, N)  # [x, z]


def assemble_full_grid(slab_dir, N, dtype, pattern):
    """Full complex grid [x, y, z] from slabs y = 0..N/2, using w(-k) = conj(w(k)) for y > N/2."""
    grid = np.zeros((N, N, N), dtype=np.complex128)
    for y in range(N // 2 + 1):
        grid[:, y, :] = read_slab(os.path.join(slab_dir, pattern.format(y=y)), N, dtype)
    idx = np.arange(N)
    mirror = (-idx) % N
    ys = np.arange(N // 2 + 1, N)
    grid[:, ys, :] = np.conj(grid[np.ix_(mirror, (-ys) % N, mirror)])
    return grid


def int_wavenumbers(N):
    k1 = np.fft.fftfreq(N) * N
    return np.meshgrid(k1, k1, k1, indexing="ij")


def test1(log, w, N, kk):
    log("=" * 72)
    log(f"TEST 1: white-noise statistics (N0 = {N}, KK limit = {kk})")
    KX, KY, KZ = int_wavenumbers(N)
    k2 = KX**2 + KY**2 + KZ**2
    nyq = (np.abs(KX) == N // 2) | (np.abs(KY) == N // 2) | (np.abs(KZ) == N // 2)
    dc = k2 == 0
    fixed = (k2 <= kk) & ~dc & ~nyq
    free = ~fixed & ~dc & ~nyq
    a2 = np.abs(w) ** 2
    log(f"  modes: free = {free.sum()}, fixed = {fixed.sum()}, Nyquist = {nyq.sum()}")
    log(f"  <|w|^2> (free modes)          = {a2[free].mean():.6f}   (expect 1)")
    log(f"  var(Re w), var(Im w) (free)   = {w.real[free].var():.6f}, {w.imag[free].var():.6f}   (expect 0.5, 0.5)")
    log(f"  <Re w>, <Im w> (free)         = {w.real[free].mean():.3e}, {w.imag[free].mean():.3e}   (expect ~0)")
    if fixed.any():
        dev = np.abs(np.sqrt(a2[fixed]) - 1.0)
        log(f"  max ||w|-1| (fixed modes)     = {dev.max():.3e}   (expect ~0)")
    log(f"  max |w| on Nyquist planes     = {np.sqrt(a2[nyq].max()) if nyq.any() else 0.0:.3e}   (expect 0)")
    log(f"  |w(k=0)|                      = {np.abs(w[0, 0, 0]):.3e}")
    kint = np.sqrt(k2)
    centers, mean_a2, counts = shell_average(a2, np.where(free, kint, -1), 16, N / 2)
    log(f"  {'k [int]':>9} {'<|w|^2>':>10} {'N_modes':>9}   (a drop near Nyquist would show the mode_weightings effect)")
    for c, a, n in zip(centers, mean_a2, counts):
        if n > 0:
            log(f"  {c:9.2f} {a:10.5f} {n:9d}")
    idx = np.arange(N)
    m = (-idx) % N
    plane = w[:, 0, :]
    herm = np.abs(plane - np.conj(plane[np.ix_(m, m)])).max()
    log(f"  max Hermitian violation, y=0  = {herm:.3e}   (expect ~0)")


def shell_average(values, kmag, nbins, kmax):
    edges = np.linspace(0.5, kmax, nbins + 1)
    which = np.digitize(kmag.ravel(), edges) - 1
    ok = (which >= 0) & (which < nbins)
    counts = np.bincount(which[ok], minlength=nbins)
    sums = np.bincount(which[ok], weights=values.ravel()[ok], minlength=nbins)
    centers = 0.5 * (edges[1:] + edges[:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        return centers, sums / counts, counts


def test2(log, dk, N, box, pk_path, kk, norm, nbins):
    log("=" * 72)
    log("TEST 2: P_ours(k) / P_input(k)")
    kin, pin = np.loadtxt(pk_path, usecols=(0, 1), unpack=True)
    KX, KY, KZ = int_wavenumbers(N)
    kint = np.sqrt(KX**2 + KY**2 + KZ**2)
    kphys = 2 * np.pi / box * kint
    valid = (kint > 0) & (np.abs(KX) < N // 2) & (np.abs(KY) < N // 2) & (np.abs(KZ) < N // 2)
    p_meas = np.abs(dk) ** 2 * norm
    p_in = np.zeros_like(kphys)
    p_in[valid] = np.exp(np.interp(np.log(kphys[valid]), np.log(kin), np.log(pin)))
    ratio = np.where(valid, p_meas / np.where(p_in > 0, p_in, 1.0), np.nan)
    k2 = kint**2
    fixed = valid & (k2 <= kk)
    if fixed.any():
        log(f"  fixed modes (k^2 <= {kk}): mean ratio = {np.nanmean(ratio[fixed]):.6f}, "
            f"max |ratio-1| = {np.nanmax(np.abs(ratio[fixed]-1)):.3e}   (expect 1, ~0)")
    centers, mean_ratio, counts = shell_average(np.nan_to_num(ratio), np.where(valid, kint, -1), nbins, N / 2)
    log(f"  {'k [int]':>9} {'k [h/Mpc]':>11} {'<P/P_in>':>10} {'N_modes':>9}")
    for c, r, n in zip(centers, mean_ratio, counts):
        if n > 0:
            log(f"  {c:9.2f} {2*np.pi/box*c:11.5f} {r:10.5f} {n:9d}")


def test3(log, dk_ours, delta_mono, N, box, nbins, kfit):
    log("=" * 72)
    log("TEST 3: cross-correlation with monofonIC delta")
    dk_mono = np.fft.fftn(delta_mono)
    KX, KY, KZ = int_wavenumbers(N)
    kint = np.sqrt(KX**2 + KY**2 + KZ**2)
    valid = (kint > 0) & (np.abs(KX) < N // 2) & (np.abs(KY) < N // 2) & (np.abs(KZ) < N // 2)
    cross = (dk_ours * np.conj(dk_mono)).real
    kv = np.where(valid, kint, -1)
    c, pxy, n = shell_average(cross, kv, nbins, N / 2)
    _, pxx, _ = shell_average(np.abs(dk_ours) ** 2, kv, nbins, N / 2)
    _, pyy, _ = shell_average(np.abs(dk_mono) ** 2, kv, nbins, N / 2)
    log(f"  {'k [int]':>9} {'k [h/Mpc]':>11} {'r(k)':>9} {'N_modes':>9}")
    for ci, a, b, d, ni in zip(c, pxy, pxx, pyy, n):
        if ni > 0:
            log(f"  {ci:9.2f} {2*np.pi/box*ci:11.5f} {a/np.sqrt(b*d):9.5f} {ni:9d}")
    # If delta_mono(x) = delta_ours(x - Delta), then arg(d_ours d_mono*) = k . Delta (Delta in grid cells).
    # Valid while |k . Delta| < pi for all fitted modes, i.e. |Delta| < N / (2 kfit).
    sel = valid & (kint <= kfit)
    phase = np.angle(dk_ours[sel] * np.conj(dk_mono[sel]))
    A = 2 * np.pi / N * np.stack([KX[sel], KY[sel], KZ[sel]], axis=1)
    delta_cells, *_ = np.linalg.lstsq(A, phase, rcond=None)
    log(f"  best-fit translation of mono relative to ours, |k| <= {kfit}: "
        f"({delta_cells[0]:+.4f}, {delta_cells[1]:+.4f}, {delta_cells[2]:+.4f}) grid cells   (expect ~0)")
    log(f"  rms residual phase = {np.sqrt(np.mean((phase - A @ delta_cells)**2)):.3e} rad")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slab-dir", help="directory with output_k_space_slab_y.<y>")
    ap.add_argument("--slab-pattern", default="output_k_space_slab_y.{y}")
    ap.add_argument("--N0", type=int, help="grid size of the slab files")
    ap.add_argument("--precision", choices=["double", "float"], default="double")
    ap.add_argument("--kk", type=float, default=1025, help="Panphasia KK limit (fixed modes k^2 <= KK)")
    ap.add_argument("--ours-slabs", help="directory with our D slabs (same layout)")
    ap.add_argument("--ours-pattern", default="D_slab_y.{y}")
    ap.add_argument("--N", type=int, help="our grid size")
    ap.add_argument("--ours-norm", type=float, default=None,
                    help="P_phys = |D|^2 * norm; default V, i.e. delta(x) = sum_k D(k) e^{ikx} "
                         "(unnormalized backward FFT of D gives delta). Use V/N^6 if D = fftn(delta).")
    ap.add_argument("--box", type=float, help="box size [Mpc/h]")
    ap.add_argument("--pk", help="input P(k) at z_start: columns k [h/Mpc], P [(Mpc/h)^3]")
    ap.add_argument("--mono-hdf5", help="monofonIC [testing] test = potentials_and_densities output")
    ap.add_argument("--mono-dset", default="delta")
    ap.add_argument("--nbins", type=int, default=16)
    ap.add_argument("--kfit", type=float, default=8, help="max integer |k| for the phase-offset fit")
    ap.add_argument("--out", default="validation_results.txt")
    args = ap.parse_args()

    log = Tee(args.out)
    log("External white-noise validation")
    log("command: " + " ".join(sys.argv))
    dtype = np.float64 if args.precision == "double" else np.float32

    if args.slab_dir:
        w = assemble_full_grid(args.slab_dir, args.N0, dtype, args.slab_pattern)
        test1(log, w, args.N0, args.kk)
        del w

    dk_ours = None
    if args.ours_slabs:
        dk_ours = assemble_full_grid(args.ours_slabs, args.N, np.float64, args.ours_pattern)
        if args.pk:
            norm = args.ours_norm if args.ours_norm is not None else args.box**3
            test2(log, dk_ours, args.N, args.box, args.pk, args.kk, norm, args.nbins)

    if args.mono_hdf5 and dk_ours is not None:
        import h5py
        with h5py.File(args.mono_hdf5, "r") as f:
            d = np.asarray(f[args.mono_dset], dtype=np.float64)
        d = d[: args.N, : args.N, : args.N]  # drop FFTW padding if present
        test3(log, dk_ours, d, args.N, args.box, args.nbins, args.kfit)

    log("=" * 72)
    log(f"results written to {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
