#!/usr/bin/env python3
"""Validation tests for loading external (monofonIC/Panphasia) white noise.

Test 1  white-noise statistics of the slab files (normalization, fixed modes, Nyquist, Hermitian y=0 plane)
Test 1b crop mapping: D(k) is a positive real multiple of w_N0 at the same signed k (needs --ours-slabs)
Test 2  P_ours(k) / P_input(k) from our linear density field (needs --ours-slabs and --pk)
Test 3  cross-correlation r(k) and phase offset between our delta and monofonIC's delta (needs --mono-hdf5)
Test 4  particle ICs (RVZel ic_* + dens_*): lattice coverage, vel/displ, displacement axes, curl,
        -i k.psi vs FFT(dens), and both vs the dumped D per mode (needs --ic-dir; --ours-slabs for the D tie)

All numbers are printed and written to --out (default: validation_results.txt).

Slab format (Willem's output_k_space_slabs branch and, for our dump, the same layout):
  one file per k_y = 0..N/2, raw complex values (re, im) in [x][z] order with z fastest, no header.

Example:
  python validate_external_noise.py --slab-dir slabs/ --N0 128 --kk 1025 \
      --ours-slabs ours/ --N 128 --box 673.2 --pk class_pk_zstart.dat \
      --mono-hdf5 output.hdf5 --ic-dir zmpi_out --out validation_results.txt
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
    k1 = np.rint(np.fft.fftfreq(N) * N)   # fftfreq * N is off by ulps (113.99999...)
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


def test_crop(log, dk_ours, w, N, N0, box, pk_path, norm):
    """Per-mode check of the k-space crop N0 -> N: D(k) = sqrt(P(k)/norm) * w_N0(k) for every
    retained mode, with w_N0 taken at the same signed integer k (index k + N0 for k < 0)."""
    log("=" * 72)
    log(f"TEST 1b: crop mapping, D(k) vs w_N0(k) at the same signed k (N0 = {N0} -> N = {N})")
    k1 = np.rint(np.fft.fftfreq(N) * N).astype(np.int64)
    src = np.where(k1 >= 0, k1, k1 + N0)
    wm = w[np.ix_(src, src, src)]
    KX, KY, KZ = int_wavenumbers(N)
    kint = np.sqrt(KX**2 + KY**2 + KZ**2)
    valid = (kint > 0) & (np.abs(KX) < N // 2) & (np.abs(KY) < N // 2) & (np.abs(KZ) < N // 2) & (wm != 0)
    ratio = dk_ours[valid] / wm[valid]
    phase = np.abs(np.angle(ratio))
    log(f"  retained modes: {valid.sum()}")
    log(f"  max |arg(D / w_N0)|           = {phase.max():.3e} rad   (expect ~1e-7: D is float32)")
    log(f"  fraction with |arg| < 1e-6    = {np.mean(phase < 1e-6):.7f}   (expect 1)")
    nyq_ours = (np.abs(KX) == N // 2) | (np.abs(KY) == N // 2) | (np.abs(KZ) == N // 2)
    log(f"  max |D| on our Nyquist planes = {np.abs(dk_ours[nyq_ours]).max():.3e}   (expect 0)")
    if pk_path:
        kin, pin = np.loadtxt(pk_path, usecols=(0, 1), unpack=True)
        kphys = 2 * np.pi / box * kint[valid]
        p_in = np.exp(np.interp(np.log(kphys), np.log(kin), np.log(pin)))
        amp = np.abs(ratio) ** 2 * norm / p_in
        log(f"  |D / w_N0|^2 * norm / P_in:   median {np.median(amp):.6f}, max |x-1| {np.max(np.abs(amp - 1)):.3e}"
            f"   (expect 1; deviations = spline vs log-linear interpolation of P)")
    # sensitivity baseline: the naive first-N-indices crop puts the wrong modes at negative k
    naive = w[:N, :N, :N][valid]
    nz = naive != 0
    naive_ok = np.mean(np.abs(np.angle(dk_ours[valid][nz] / naive[nz])) < 1e-6)
    log(f"  same phase test with the naive w[:N,:N,:N] crop: fraction {naive_ok:.4f}"
        f"   (1 when N0 = N; ~1/8 when N0 > N, only all-positive k match)")


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


RVZEL_DTYPE = np.dtype([("i", "<u2"), ("j", "<u2"), ("k", "<u2"), ("pad", "<u2"),
                        ("displ", "<f4", 3), ("vel", "<f4", 3)])


def read_rvzel_lattice(ic_dir, N, ic_glob):
    """Displacements [3, N, N, N] and dens [N, N, N] on the lattice indexed by the stored (i, j, k).

    dens_<n> (float32, one value per particle, same record order as ic_<n>) is optional.
    Velocities are not stored; their proportionality to the displacements is reduced on the fly.
    """
    import glob
    files = sorted(f for f in glob.glob(os.path.join(ic_dir, ic_glob)) if os.path.basename(f)[3:].isdigit())
    if not files:
        raise FileNotFoundError(f"no {ic_glob} files in {ic_dir}")
    disp = np.zeros((3, N, N, N), dtype=np.float64)
    dens = np.zeros((N, N, N), dtype=np.float64)
    hits = np.zeros((N, N, N), dtype=np.uint8)
    have_dens = True
    sum_vd = sum_dd = 0.0
    npart = npad = 0
    vel_max_dev = 0.0
    for f in files:
        p = np.fromfile(f, dtype=RVZEL_DTYPE)
        i, j, k = (p[c].astype(np.intp) for c in ("i", "j", "k"))
        if max(i.max(), j.max(), k.max()) >= N:
            raise ValueError(f"{f}: lattice index >= N = {N}")
        d = p["displ"].astype(np.float64)
        v = p["vel"].astype(np.float64)
        disp[:, i, j, k] = d.T
        np.add.at(hits, (i, j, k), 1)
        sum_vd += float(np.sum(v * d))
        sum_dd += float(np.sum(d * d))
        npart += p.size
        npad += int(np.count_nonzero(p["pad"]))
        fd = os.path.join(os.path.dirname(f), "dens_" + os.path.basename(f)[3:])
        if have_dens and os.path.exists(fd):
            dv = np.fromfile(fd, dtype="<f4")
            if dv.size != p.size:
                raise ValueError(f"{fd}: {dv.size} values for {p.size} particles in {f}")
            dens[i, j, k] = dv
        else:
            have_dens = False
    vfac = sum_vd / sum_dd if sum_dd > 0 else np.nan
    for f in files:
        p = np.fromfile(f, dtype=RVZEL_DTYPE)
        d = p["displ"].astype(np.float64)
        vel_max_dev = max(vel_max_dev, float(np.abs(p["vel"] - vfac * d).max()))
    info = dict(nfiles=len(files), npart=npart, npad=npad, vfac=vfac, vel_max_dev=vel_max_dev,
                disp_rms=float(np.sqrt(disp.var(axis=(1, 2, 3)).mean())))
    return disp, (dens if have_dens else None), hits, info


def complex_corr(a, b):
    return float(np.real(np.vdot(b, a)) / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real))


def fit_scale(a, b):
    """Complex c minimizing |a - c b|^2."""
    return np.vdot(b, a) / np.vdot(b, b).real


def log_rel_residual(log, label, resid, ref_abs, kint, valid, nbins, N, floor=None):
    rel2 = np.where(valid, (resid / np.where(ref_abs > 0, ref_abs, 1.0)) ** 2, 0.0)
    centers, mean_rel2, counts = shell_average(rel2, np.where(valid, kint, -1), nbins, N / 2)
    if floor is not None:
        _, mean_floor2, _ = shell_average(np.where(valid, floor, 0.0), np.where(valid, kint, -1), nbins, N / 2)
    log(f"  {label}: max per-mode relative residual = {np.sqrt(rel2[valid].max()):.3e} "
        f"(all non-DC, non-Nyquist modes, incl. corners)")
    head = f"  {'k [int]':>9} {'rms rel':>11}" + (f" {'f32 floor':>11}" if floor is not None else "") + f" {'N_modes':>9}"
    log(head)
    for idx, (c, r2, n) in enumerate(zip(centers, mean_rel2, counts)):
        if n > 0:
            line = f"  {c:9.2f} {np.sqrt(r2):11.3e}"
            if floor is not None:
                line += f" {np.sqrt(mean_floor2[idx]):11.3e}"
            log(line + f" {n:9d}")


def float32_floor(field_k):
    """|FFT(float32(real(IFFT(field_k)))) - field_k|^2 / |field_k|^2: the float32 storage floor per mode."""
    x = np.fft.ifftn(field_k).real
    back = np.fft.fftn(x.astype(np.float32).astype(np.float64))
    return np.abs(back - field_k) ** 2 / np.maximum(np.abs(field_k) ** 2, 1e-300)


def test4(log, dk_ours, ic_dir, ic_glob, N, box, nbins):
    """End-to-end check of the particle ICs (RVZel ic_* + dens_*) against the dumped D."""
    from itertools import permutations
    log("=" * 72)
    log(f"TEST 4: particle ICs in {ic_dir} (RVZel, {RVZEL_DTYPE.itemsize} bytes/particle)")
    disp, dens, hits, info = read_rvzel_lattice(ic_dir, N, ic_glob)
    missing = int(np.count_nonzero(hits == 0))
    dup = int(np.count_nonzero(hits > 1))
    log(f"  files = {info['nfiles']}, particles = {info['npart']} (expect {N**3}), "
        f"lattice sites missing = {missing}, duplicated = {dup}   (expect 0, 0)")
    log(f"  non-zero padding words = {info['npad']}   (expect 0)")
    log(f"  vel = c * displ: fitted c = {info['vfac']:.7f}, max |vel - c displ| = {info['vel_max_dev']:.3e}   "
        f"(c = (sqrt(1+24 f_cluster)-1)/4, = 1 for f_cluster = 1)")
    log(f"  rms displacement per component = {info['disp_rms']:.4e} (file units)")

    KX, KY, KZ = int_wavenumbers(N)
    kvec = (KX, KY, KZ)
    k2 = KX**2 + KY**2 + KZ**2
    kint = np.sqrt(k2)
    valid = (k2 > 0) & (np.abs(KX) < N // 2) & (np.abs(KY) < N // 2) & (np.abs(KZ) < N // 2)
    kphys = 2 * np.pi / box
    psi = [np.fft.fftn(disp[c]) for c in range(3)]
    del disp

    log("  -- (a) particle frame: axes = stored (i, j, k), which Abacus maps to (x, y, z)")
    if dens is not None:
        delta = np.fft.fftn(dens)
    else:
        delta = sum(-1j * kphys * kvec[c] * psi[c] for c in range(3))
        log("  no dens_* files: using -i k.psi as the density reference")
    ksafe = np.where(k2 > 0, k2, 1.0)
    log("  corr(psi_c, i k_a delta / k^2), rows c = displ[0..2], columns a = axis i, j, k (expect diagonal +-1):")
    for c in range(3):
        row = [complex_corr(psi[c][valid], (1j * kvec[a] / ksafe * delta)[valid]) for a in range(3)]
        log("    displ[{}]  ".format(c) + "  ".join(f"{r:+9.6f}" for r in row))
    kpsi = sum(kvec[c] * psi[c] for c in range(3))
    kxpsi2 = (np.abs(KY * psi[2] - KZ * psi[1]) ** 2 + np.abs(KZ * psi[0] - KX * psi[2]) ** 2
              + np.abs(KX * psi[1] - KY * psi[0]) ** 2)
    curl_frac = np.sqrt(kxpsi2[valid].sum() / (np.abs(kpsi[valid]) ** 2).sum())
    log(f"  irrotational: sqrt(sum |k x psi|^2 / sum |k . psi|^2) = {curl_frac:.3e}   (expect ~float32 noise)")
    theta = -1j * kphys * kpsi
    if dens is not None:
        b = fit_scale(theta[valid], delta[valid])
        log(f"  -i k.psi = b * FFT(dens): b = {b.real:+.7f} {b.imag:+.2e}i   "
            f"(expect +1 if displ is in Mpc/h and dens = -div psi)")
        log_rel_residual(log, "-i k.psi vs b FFT(dens)", np.abs(theta - b * delta), np.abs(b * delta),
                         kint, valid, nbins, N)

    if dk_ours is None:
        log("  -- (b) skipped: no --ours-slabs")
        return
    log("  -- (b) tie to the dumped D (frame [x, y, z] of the D slabs)")
    best = None
    for perm in permutations(range(3)):
        dt = np.transpose(dk_ours, perm)
        for conj in (False, True):
            ref = np.conj(dt) if conj else dt
            r = complex_corr(delta[valid], ref[valid])
            if best is None or abs(r) > abs(best[0]):
                best = (r, perm, conj)
    r, perm, conj = best
    names = "xyz"
    log(f"  best match: particle axes (i, j, k) = D axes ({', '.join(names[p] for p in perm)}), "
        f"{'conj(D)' if conj else 'D'}, r = {r:+.8f}")
    ref = np.transpose(dk_ours, perm)
    ref = np.conj(ref) if conj else ref
    a = fit_scale(delta[valid], ref[valid])
    log(f"  FFT(dens) = A * D: A / N^3 = {a.real / N**3:+.7f} {a.imag / N**3:+.2e}i   "
        f"(expect +1: dens = unnormalized backward FFT of D)")
    floor = float32_floor(a * ref)
    log_rel_residual(log, "FFT(dens) vs A D", np.abs(delta - a * ref), np.abs(a * ref),
                     kint, valid, nbins, N, floor=floor)
    at = fit_scale(theta[valid], ref[valid])
    log(f"  -i k.psi = A' * D: A' / N^3 = {at.real / N**3:+.7f} {at.imag / N**3:+.2e}i")
    log_rel_residual(log, "-i k.psi vs A' D", np.abs(theta - at * ref), np.abs(at * ref),
                     kint, valid, nbins, N)


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
    ap.add_argument("--mono-sign", type=float, choices=[1.0, -1.0], default=-1.0,
                    help="factor applied to the monofonIC field; the testing-mode 'delta' is laplace(phi), "
                         "the opposite sign to our delta, so the default -1 makes r(k) = +1")
    ap.add_argument("--ic-dir", help="Zeldovich-MPI particle output (RVZel ic_* and optional dens_*) for test 4")
    ap.add_argument("--ic-glob", default="ic_*")
    ap.add_argument("--nbins", type=int, default=16)
    ap.add_argument("--kfit", type=float, default=8, help="max integer |k| for the phase-offset fit")
    ap.add_argument("--out", default="validation_results.txt")
    args = ap.parse_args()

    log = Tee(args.out)
    log("External white-noise validation")
    log("command: " + " ".join(sys.argv))
    dtype = np.float64 if args.precision == "double" else np.float32

    w = None
    if args.slab_dir:
        w = assemble_full_grid(args.slab_dir, args.N0, dtype, args.slab_pattern)
        test1(log, w, args.N0, args.kk)

    dk_ours = None
    if args.ours_slabs:
        dk_ours = assemble_full_grid(args.ours_slabs, args.N, np.float64, args.ours_pattern)
        norm = args.ours_norm if args.ours_norm is not None else args.box**3
        if w is not None:
            test_crop(log, dk_ours, w, args.N, args.N0, args.box, args.pk, norm)
        if args.pk:
            test2(log, dk_ours, args.N, args.box, args.pk, args.kk, norm, args.nbins)
    del w

    if args.mono_hdf5 and dk_ours is not None:
        import h5py
        with h5py.File(args.mono_hdf5, "r") as f:
            d = np.asarray(f[args.mono_dset], dtype=np.float64)
        d = args.mono_sign * d[: args.N, : args.N, : args.N]  # drop FFTW padding if present
        log(f"monofonIC field multiplied by --mono-sign = {args.mono_sign:+.0f}")
        test3(log, dk_ours, d, args.N, args.box, args.nbins, args.kfit)

    if args.ic_dir:
        test4(log, dk_ours, args.ic_dir, args.ic_glob, args.N, args.box, args.nbins)

    log("=" * 72)
    log(f"results written to {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
