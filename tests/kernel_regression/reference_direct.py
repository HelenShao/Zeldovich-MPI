#!/usr/bin/env python3
"""Independent single-process reference for the f_NL kernels (small N only).

Re-implements the kernel math directly with numpy on the full field: no pyfftw, no
rank/batch/pencil emulation. It is a second reference next to kernel_MPI.py, so a bug
shared by both the C code and kernel_MPI.py's MPI emulation would still show up.

Conventions follow kernel_MPI.py exactly: Hermitian completion of the y_primary files,
zero padding to M (frequency n mod N -> n mod M), z = (G + iH) phi packed transform,
Q(x) = Im(z^2)/2 / M^3 with unnormalized transforms, Q(k=0) = 0, crop back to N, times J.

Writes phi_output_rank{r}_y{y}.bin with the same ownership as the C code.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_test_phi import rank_pair_range  # noqa: E402


def load_phi(bin_dir, N, num_ranks):
    phi = np.zeros((N, N, N), dtype=np.complex128)   # [x, y, z]
    for r in range(num_ranks):
        for y in rank_pair_range(r, N // 2 + 1, num_ranks):
            s = np.fromfile(os.path.join(bin_dir, f"Phi_gaussian_rank{r}_y{y}.bin"), dtype=np.complex128)
            phi[:, y, :] = s.reshape(N, N).T
    Nh = N // 2
    neg = (-np.arange(N)) % N
    # y mirrors
    for y in range(1, Nh):
        phi[:, N - y, :] = np.conj(phi[np.ix_(neg, [y], neg)][:, 0, :])
    # in-plane fix-ups of the self-conjugate planes (same order as kernel_MPI.py)
    for y0 in (0, Nh):
        P = phi[:, y0, :]
        zi = np.arange(1, Nh)
        P[np.ix_(neg, N - zi)] = np.conj(P[:, zi])
        xi = np.arange(1, Nh)
        P[N - xi, 0] = np.conj(P[xi, 0])
        P[N - xi, Nh] = np.conj(P[xi, Nh])
        for xx in (0, Nh):
            for zz in (0, Nh):
                P[xx, zz] = P[xx, zz].real
    return phi


def pad_index(N, M):
    i = np.arange(N)
    return np.where(i < N // 2, i, i + (M - N))


def kpow(k, p):
    if p == 0:
        return np.ones_like(k)
    out = np.zeros_like(k)
    nz = k != 0
    out[nz] = k[nz] ** p
    return out


SHAPES = {  # name: (pG, pH, pJ, cJ)
    "A1": (0, 0, 0, 1.0), "C1": (1, 1, -2, 1 / 3), "B1": (1, 0, -1, 1.0), "B2": (2, 0, -2, 1.0),
    "D1": (-1, -1, 2, 1.0), "E2": (2, -1, -1, 1 / 2), "F1": (-1, 0, 1, 1.0), "G1": (-1, 1, 0, 1.0),
}


def terms_for(fnl):
    if fnl == "local":
        return [("A1", 1.0)]
    if fnl == "equil":
        return [("A1", -3.0), ("C1", -6.0), ("B1", 4.0), ("B2", 2.0)]
    if fnl == "orthoCMB":
        return [("A1", -9.0), ("C1", -24.0), ("B1", 10.0), ("B2", 8.0)]
    p = 27 / (-21 + (743) / (7 * (20 * np.pi**2 - 193)))
    t = (2 + (60 * p) / 27) / (6 * (1 + (15 * p) / 27))
    return [("A1", -(1 + 9 * p / 27)), ("C1", -(2 + 60 * p / 27)),
            ("B1", (1 + 15 * p / 27) * (1 - t)), ("B2", (1 + 15 * p / 27) * t),
            ("D1", p / 27), ("E2", -(20 * p) / 27), ("F1", 6 * p / 27), ("G1", 15 * p / 27)]


def kernel_output(phi, L, M, fnl):
    N = phi.shape[0]
    n = np.fft.fftfreq(N) * N
    kx, ky, kz = np.meshgrid(n, n, n, indexing="ij")
    k = (2 * np.pi / L) * np.sqrt(kx**2 + ky**2 + kz**2)
    pi = pad_index(N, M)
    out = np.zeros_like(phi)
    for name, coef in terms_for(fnl):
        pG, pH, pJ, cJ = SHAPES[name]
        zk = np.zeros((M, M, M), dtype=np.complex128)
        zk[np.ix_(pi, pi, pi)] = kpow(k, pG) * phi + 1j * (kpow(k, pH) * phi)
        zx = np.fft.ifftn(zk) * M**3             # unnormalized backward transform
        q = zx.real * zx.imag / M**3             # = Im(z^2)/2 / M^3
        Qk = np.fft.fftn(q)
        Qk[0, 0, 0] = 0
        out += coef * (Qk[np.ix_(pi, pi, pi)] * (cJ * kpow(k, pJ)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ranks", type=int, required=True)
    ap.add_argument("--box-size", type=float, required=True)
    ap.add_argument("--fnl", required=True, choices=["local", "equil", "orthoCMB", "orthoLSS"])
    ap.add_argument("--M", type=int, help="padded size (default 2N)")
    ap.add_argument("--zero-nyquist", action="store_true", help="zero all kx/ky/kz = N/2 modes (KERNEL_ZERO_NYQUIST=1)")
    args = ap.parse_args()

    s0 = np.fromfile(os.path.join(args.bin_dir, "Phi_gaussian_rank0_y0.bin"), dtype=np.complex128)
    N = int(round(np.sqrt(s0.size)))
    M = args.M or 2 * N
    phi = load_phi(args.bin_dir, N, args.ranks)
    if args.zero_nyquist:
        phi[N // 2, :, :] = 0
        phi[:, N // 2, :] = 0
        phi[:, :, N // 2] = 0
    out = kernel_output(phi, args.box_size, M, args.fnl)
    os.makedirs(args.out_dir, exist_ok=True)
    for r in range(args.ranks):
        for y in rank_pair_range(r, N // 2 + 1, args.ranks):
            np.ascontiguousarray(out[:, y, :].T).tofile(os.path.join(args.out_dir, f"phi_output_rank{r}_y{y}.bin"))


if __name__ == "__main__":
    main()
