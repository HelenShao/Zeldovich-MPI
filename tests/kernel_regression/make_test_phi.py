#!/usr/bin/env python3
"""Write synthetic Phi_gaussian_y.{y} inputs (+ a minimal par2) for Kernel_MPI tests.

The field is a seeded real Gaussian random field shaped by a power law, so its FFT is
exactly Hermitian, like the Zeldovich_MPI dumps. Slices y in [0, N/2] are written keyed
by y only, the same names ZD_MPI_generation.c dumps into ZD_PNG_DIR.
File layout: N*N complex128, x fastest then z (element (x, z) at z*N + x).

The field depends only on (N, seed, slope), so one input directory serves every rank count.
"""
import argparse
import os

import numpy as np

FNL_KEYS = {
    "local": "ZD_f_NL_local",
    "equil": "ZD_f_NL_equil",
    "orthoCMB": "ZD_f_NL_orthoCMB",
    "orthoLSS": "ZD_f_NL_orthoLSS",
}


def make_phi_k(N, seed, slope=-1.5):
    """Hermitian phi(k) on an N^3 grid, indexed [x, y, z]."""
    rng = np.random.default_rng(seed)
    white = rng.standard_normal((N, N, N))
    n = np.fft.fftfreq(N) * N
    kx, ky, kz = np.meshgrid(n, n, n, indexing="ij")
    kmag = np.sqrt(kx**2 + ky**2 + kz**2)
    kmag[0, 0, 0] = 1.0
    amp = kmag**slope
    amp[0, 0, 0] = 0.0
    return np.fft.fftn(white) * amp / N**1.5


def write_par2(path, box_size, num_z_ranks, fnl, png_dir=None):
    with open(path, "w") as f:
        f.write("# Kernel_MPI regression test parameters\n")
        f.write(f"BoxSize = {box_size!r}\n")
        f.write(f"NumZRanks = {num_z_ranks}\n")
        if png_dir is not None:
            f.write(f'ZD_PNG_DIR = "{png_dir}"\n')
        for name, key in FNL_KEYS.items():
            f.write(f"{key} = {1.0 if name == fnl else 0.0}\n")


def write_inputs(out_dir, N, seed, slope=-1.5):
    os.makedirs(out_dir, exist_ok=True)
    phi = make_phi_k(N, seed, slope)
    for y in range(N // 2 + 1):
        plane_zx = np.ascontiguousarray(phi[:, y, :].T, dtype=np.complex128)
        plane_zx.tofile(os.path.join(out_dir, f"Phi_gaussian_y.{y}"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--num-z-ranks", type=int, default=1)
    ap.add_argument("--fnl", choices=list(FNL_KEYS), default="local")
    ap.add_argument("--box-size", type=float, default=500.0)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--slope", type=float, default=-1.5)
    ap.add_argument("--out", required=True, help="directory for the Phi_gaussian_y.{y} slices")
    ap.add_argument("--par", help="par2 path to write (default: <out>/param.par2)")
    args = ap.parse_args()

    assert args.N % 2 == 0
    write_inputs(args.out, args.N, args.seed, args.slope)
    write_par2(args.par or os.path.join(args.out, "param.par2"), args.box_size, args.num_z_ranks, args.fnl,
               png_dir=os.path.abspath(args.out))


if __name__ == "__main__":
    main()
