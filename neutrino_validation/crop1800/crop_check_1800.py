#!/usr/bin/env python3
"""Crop test (b) at N = 1800 from the FLAMINGO noise (N0 = 2048), flat P(k) run.

1. k-space, per y-plane: D_slab_y.<y> (N = 1800) vs the N0 = 2048 noise slab of the same k_y,
   taken at the same signed (k_x, k_z). With P = 1/L^3, D = w / L^1.5 exactly (up to float32).
2. real space: a sub-block of our dens_* (the first NSUB lattice planes along the ic_ file axis)
   vs the same planes of monofonIC's 1800^3 white_noise field (its own crop of the N0 grid).
   Axis order and sign conventions differ (Abacus (x, y, z) = monofonIC (z, y, x)), so all
   candidate orientations are correlated and the best is reported; expect |r| ~ 1 for one of
   them (slightly below 1: our Nyquist planes are zeroed) and ~0 for the others.
"""
import argparse
import glob
import os
import sys

import numpy as np

RVZEL = np.dtype([("i", "<u2"), ("j", "<u2"), ("k", "<u2"), ("pad", "<u2"),
                  ("displ", "<f4", 3), ("vel", "<f4", 3)])


def read_plane(path, n, dtype=np.float64):
    raw = np.fromfile(path, dtype=dtype)
    if raw.size != 2 * n * n:
        raise ValueError(f"{path}: {raw.size} values, expected {2 * n * n}")
    return (raw[0::2] + 1j * raw[1::2]).reshape(n, n)  # [x, z]


def kspace_check(log, d_dir, noise_dir, N, N0, box, ys):
    k1 = np.rint(np.fft.fftfreq(N) * N).astype(np.int64)   # plain astype truncates 113.99999... to 113
    src = np.where(k1 >= 0, k1, k1 + N0)
    KX, KZ = np.meshgrid(k1, k1, indexing="ij")
    expect = box ** -1.5
    ok = True
    log(f"[k-space] D / w_N0 per plane, expect the real constant L^-1.5 = {expect:.6e}")
    for y in ys:
        D = read_plane(os.path.join(d_dir, f"D_slab_y.{y}"), N)
        w = read_plane(os.path.join(noise_dir, f"output_k_space_slab_y.{y}"), N0)[np.ix_(src, src)]
        ky = y if y <= N // 2 else y - N
        valid = (np.abs(KX) < N // 2) & (np.abs(KZ) < N // 2) & (abs(ky) < N // 2) & (w != 0)
        if y == 0:
            valid &= ~((KX == 0) & (KZ == 0))
        r = D[valid] / w[valid]
        rel = np.max(np.abs(r / expect - 1)) if r.size else 0.0
        nyq = ~((np.abs(KX) < N // 2) & (np.abs(KZ) < N // 2)) | (abs(ky) >= N // 2)
        dnyq = np.abs(D[nyq]).max() if nyq.any() else 0.0
        good = r.size > 0 and rel < 2e-7 and dnyq == 0 if abs(ky) < N // 2 else dnyq == 0
        ok &= bool(good)
        log(f"  y={y:4d}: modes {r.size:8d}  max|D/(w L^-1.5) - 1| = {rel:.2e}  max|D| on Nyquist = {dnyq:.1e}"
            f"  -> {'PASS' if good else 'FAIL'}")
    return ok


def noise_stats(log, noise_dir, N0):
    """Every slab of the N0 noise: size, <|w|^2> per plane (expect ~1, no empty or garbage planes),
    Hermitian symmetry of the k_y = 0 plane."""
    ny = N0 // 2 + 1
    means = np.zeros(ny)
    for y in range(ny):
        w = read_plane(os.path.join(noise_dir, f"output_k_space_slab_y.{y}"), N0)
        nz = w != 0
        means[y] = np.mean(np.abs(w[nz]) ** 2) if nz.any() else 0.0
        if y == 0:
            i = (-np.arange(N0)) % N0
            herm = np.max(np.abs(w - np.conj(w[np.ix_(i, i)])))
    # the Nyquist plane y = N0/2 is zeroed by Panphasia's mode weightings; exclude it from the range
    body = means[:-1]
    lo, hi = body.min(), body.max()
    ok = bool(0.9 < lo and hi < 1.1 and herm < 1e-12)
    log(f"[noise] {ny} slabs of N0 = {N0}: per-plane <|w|^2> in [{lo:.4f}, {hi:.4f}] "
        f"(planes 0..{ny - 2}), mean {body.mean():.5f}, y = {ny - 1}: {means[-1]:.3g}; "
        f"Hermitian violation y=0: {herm:.2e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def lowk_check(log, noise_dir, N0, ref_dir, N0_ref, kmax=128, kpass=32):
    """Panphasia is multi-resolution: the same descriptor at a lower N0 must give the same low-k modes
    (both slab sets have unit variance per mode). Correlate the N0 slabs with a lower-N0 reference
    at the same signed k, |k_x|, |k_y|, |k_z| <= kmax. Pass: r > 0.99 in every shell with k <= kpass."""
    k1 = np.arange(-kmax, kmax + 1)
    ia, ib = k1 % N0, k1 % N0_ref
    nb = int(np.ceil(np.sqrt(3) * kmax))
    sab, saa, sbb, cnt = (np.zeros(nb + 1, dtype=c) for c in (np.complex128, float, float, np.int64))
    KX, KZ = np.meshgrid(k1, k1, indexing="ij")
    for y in range(kmax + 1):
        a = read_plane(os.path.join(noise_dir, f"output_k_space_slab_y.{y}"), N0)[np.ix_(ia, ia)]
        b = read_plane(os.path.join(ref_dir, f"output_k_space_slab_y.{y}"), N0_ref)[np.ix_(ib, ib)]
        kk = np.sqrt(KX ** 2 + y ** 2 + KZ ** 2)
        sel = (kk > 0) & (kk <= kmax)
        idx = np.rint(kk[sel]).astype(np.int64)
        np.add.at(sab, idx, a[sel] * np.conj(b[sel]))
        np.add.at(saa, idx, np.abs(a[sel]) ** 2)
        np.add.at(sbb, idx, np.abs(b[sel]) ** 2)
        np.add.at(cnt, idx, 1)
    log(f"[low k] N0 = {N0} slabs vs the same descriptor at N0 = {N0_ref} ({ref_dir})")
    log("    k [int]   N_modes   Re r(k)   Im r(k)   <|w|^2> ratio")
    ok = True
    for k in [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]:
        if k > kmax or cnt[k] == 0:
            continue
        r = sab[k] / np.sqrt(saa[k] * sbb[k])
        good = r.real > 0.99 if k <= kpass else True
        ok &= bool(good)
        log(f"    {k:7d} {cnt[k]:9d}  {r.real:8.5f}  {r.imag:+8.1e}  {saa[k] / sbb[k]:8.5f}"
            f"{'' if k <= kpass else '   (informational)'}")
    log(f"  r > 0.99 for k <= {kpass} -> {'PASS' if ok else 'FAIL'}")
    return ok


def read_ours_subblock(ic_dir, N, nsub):
    """Our dens_ values on the first nsub lattice planes along the ic_ file axis, as [nsub, N, N]."""
    files = sorted(f for f in glob.glob(os.path.join(ic_dir, "ic_*")) if os.path.basename(f)[3:].isdigit())
    I, J, K, V = [], [], [], []
    axis = None
    for f in files:
        rec = np.fromfile(f, RVZEL)
        dens = np.fromfile(os.path.join(ic_dir, "dens_" + os.path.basename(f)[3:]), np.float32)
        if rec.size == 0:
            continue
        idx = np.stack([rec["i"], rec["j"], rec["k"]]).astype(np.int64)
        if axis is None:
            axis = int(np.argmin([np.ptp(a) for a in idx]))
        sel = idx[axis] < nsub
        if not sel.any():
            break
        I.append(idx[0][sel]); J.append(idx[1][sel]); K.append(idx[2][sel]); V.append(dens[sel])
    I, J, K, V = (np.concatenate(a) for a in (I, J, K, V))
    idx = [I, J, K]
    order = [axis] + [a for a in range(3) if a != axis]
    block = np.full((nsub, N, N), np.nan, dtype=np.float64)
    block[idx[order[0]], idx[order[1]], idx[order[2]]] = V
    return block, axis


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float(np.sum(a * b) / np.sqrt(np.sum(a * a) * np.sum(b * b)))


def realspace_check(log, ic_dir, mono_h5, N, nsub):
    import h5py
    ours, axis = read_ours_subblock(ic_dir, N, nsub)
    nmiss = int(np.isnan(ours).sum())
    log(f"[real space] our dens_ block: first {nsub} planes along RVZel index {'ijk'[axis]}, missing {nmiss}")
    if nmiss:
        return False
    with h5py.File(mono_h5, "r") as f:
        dset = f["white_noise"]
        first = np.asarray(dset[:nsub, :N, :N], dtype=np.float64)          # mono axis 0 = planes
        mid = np.asarray(dset[:N, :nsub, :N], dtype=np.float64)            # mono axis 1 = planes
        last = np.asarray(dset[:N, :N, :nsub], dtype=np.float64)           # mono axis 2 = planes
    cands = {
        "mono[p,a,b]": first,
        "mono[p,b,a]": first.transpose(0, 2, 1),
        "mono[a,p,b] -> [p,a,b]": mid.transpose(1, 0, 2),
        "mono[a,p,b] -> [p,b,a]": mid.transpose(1, 2, 0),
        "mono[a,b,p] -> [p,b,a]": last.transpose(2, 1, 0),
        "mono[a,b,p] -> [p,a,b]": last.transpose(2, 0, 1),
    }
    best = 0.0
    for name, m in cands.items():
        r = corr(ours, m)
        best = max(best, abs(r))
        log(f"  r(ours, {name:24s}) = {r:+.6f}")
    ok = best > 0.99
    log(f"  best |r| = {best:.6f} (expect > 0.99; Nyquist planes zeroed in ours) -> {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--N", type=int, default=1800)
    p.add_argument("--N0", type=int, default=2048)
    p.add_argument("--box", type=float, default=673.2)
    p.add_argument("--d-dir", default="D_slabs")
    p.add_argument("--noise-dir", required=True)
    p.add_argument("--ic-dir", default="ic")
    p.add_argument("--mono-h5", required=True, help="skipped if missing (monofonIC's white_noise output)")
    p.add_argument("--ref-noise-dir", default=None, help="same descriptor at a lower N0 (low-k check)")
    p.add_argument("--ref-N0", type=int, default=512)
    p.add_argument("--nsub", type=int, default=64)
    p.add_argument("--ys", default="0,1,2,449,450,451,898,899,900")
    p.add_argument("--out", default="crop1800_results.txt")
    a = p.parse_args()

    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"crop test (b): N = {a.N} from N0 = {a.N0}, box {a.box}")
    oks = [noise_stats(log, a.noise_dir, a.N0)]
    if a.ref_noise_dir:
        oks.append(lowk_check(log, a.noise_dir, a.N0, a.ref_noise_dir, a.ref_N0))
    oks.append(kspace_check(log, a.d_dir, a.noise_dir, a.N, a.N0, a.box, [int(y) for y in a.ys.split(",")]))
    if os.path.exists(a.mono_h5):
        oks.append(realspace_check(log, a.ic_dir, a.mono_h5, a.N, a.nsub))
    else:
        log(f"[real space] skipped: {a.mono_h5} not found")
    ok = all(oks)
    log(f"OVERALL: {'PASS' if ok else 'FAIL'}")
    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
