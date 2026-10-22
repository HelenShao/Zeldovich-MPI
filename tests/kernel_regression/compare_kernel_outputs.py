#!/usr/bin/env python3
"""Compare two directories of kernel_y.{y} kernel outputs.

Slices are matched by y, so outputs from runs with different rank counts can be compared. Exit status is 0 if all tolerances pass, 1 otherwise.

Reported metrics (B is the reference):
  rel_l2   ||A - B|| / ||B|| over all slices
  max_abs  max |A - B|
  max_rel  max |A - B| / |B| over modes with |B| > 1e-8 * max|B|
  herm     Hermitian residual of the y=0 plane of A (Nyquist rows/columns excluded),
           relative to max|A|; a sanity check that the output is a real field's FFT
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

_PAT = re.compile(r"kernel_y\.(\d+)$")


def load_outputs(directory):
    """Return ({y: (N, N) complex array indexed [z, x]}, N)."""
    slices = {}
    for path in glob.glob(os.path.join(directory, "kernel_y.*")):
        m = _PAT.search(path)
        if not m:
            continue
        y = int(m.group(1))
        data = np.fromfile(path, dtype=np.complex128)
        N = int(round(np.sqrt(data.size)))
        if N * N != data.size:
            raise SystemExit(f"{path}: {data.size} values is not N*N")
        slices[y] = data.reshape(N, N)
    if not slices:
        raise SystemExit(f"no kernel_y.* files in {directory}")
    Ns = {s.shape[0] for s in slices.values()}
    if len(Ns) != 1:
        raise SystemExit(f"{directory}: slices of different sizes {sorted(Ns)}")
    return slices, Ns.pop()


def hermitian_residual(plane):
    """max |P(x,z) - conj P(-x,-z)| on the y=0 plane, excluding the Nyquist row/column."""
    N = plane.shape[0]
    idx = (-np.arange(N)) % N
    mirrored = np.conj(plane[np.ix_(idx, idx)])
    keep = np.ones(N, dtype=bool)
    keep[N // 2] = False
    diff = np.abs(plane - mirrored)[np.ix_(keep, keep)]
    scale = np.max(np.abs(plane))
    return float(np.max(diff) / scale) if scale > 0 else 0.0


def compare(dir_a, dir_b, rel_l2_tol=1e-10, max_rel_tol=1e-8):
    a, Na = load_outputs(dir_a)
    b, Nb = load_outputs(dir_b)
    result = {"dir_a": dir_a, "dir_b": dir_b, "N": Nb, "num_slices": len(b), "failures": []}
    if Na != Nb:
        result["failures"].append(f"N differs: {Na} vs {Nb}")
        return result
    if sorted(a) != sorted(b):
        result["failures"].append(f"y sets differ: only in A {sorted(set(a) - set(b))}, only in B {sorted(set(b) - set(a))}")
        return result

    num = den = 0.0
    max_abs = max_rel = 0.0
    worst = None
    gmax = max(np.max(np.abs(s)) for s in b.values())
    for y in sorted(b):
        d = np.abs(a[y] - b[y])
        num += float(np.sum(d**2))
        den += float(np.sum(np.abs(b[y]) ** 2))
        mask = np.abs(b[y]) > 1e-8 * gmax
        rel = np.where(mask, d / np.where(mask, np.abs(b[y]), 1.0), 0.0)
        if np.max(d) > max_abs:
            max_abs = float(np.max(d))
        if np.max(rel) > max_rel:
            max_rel = float(np.max(rel))
            z, x = np.unravel_index(np.argmax(rel), rel.shape)
            worst = {"y": y, "x": int(x), "z": int(z), "a": str(a[y][z, x]), "b": str(b[y][z, x])}

    result["rel_l2"] = float(np.sqrt(num / den)) if den > 0 else float(np.sqrt(num))
    result["max_abs"] = max_abs
    result["max_rel"] = max_rel
    result["worst"] = worst
    result["herm_y0"] = hermitian_residual(a[0]) if 0 in a else None

    if not result["rel_l2"] <= rel_l2_tol:
        result["failures"].append(f"rel_l2 {result['rel_l2']:.3e} > {rel_l2_tol:.1e}")
    if not result["max_rel"] <= max_rel_tol:
        result["failures"].append(f"max_rel {result['max_rel']:.3e} > {max_rel_tol:.1e} at {worst}")
    if result["herm_y0"] is not None and not result["herm_y0"] <= 1e-10:
        result["failures"].append(f"y=0 plane not Hermitian (residual {result['herm_y0']:.3e})")
    return result


def summary_line(r):
    if "rel_l2" not in r:
        return "; ".join(r["failures"])
    return (f"rel_l2={r['rel_l2']:.2e} max_rel={r['max_rel']:.2e} max_abs={r['max_abs']:.2e} "
            f"herm_y0={r['herm_y0']:.1e}" if r["herm_y0"] is not None else
            f"rel_l2={r['rel_l2']:.2e} max_rel={r['max_rel']:.2e} max_abs={r['max_abs']:.2e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir_a", help="output directory under test (e.g. the C run)")
    ap.add_argument("dir_b", help="reference output directory (e.g. kernel_MPI.py)")
    ap.add_argument("--rel-l2-tol", type=float, default=1e-10)
    ap.add_argument("--max-rel-tol", type=float, default=1e-8)
    ap.add_argument("--json", help="write the result dict here")
    args = ap.parse_args()

    r = compare(args.dir_a, args.dir_b, args.rel_l2_tol, args.max_rel_tol)
    print(("PASS  " if not r["failures"] else "FAIL  ") + summary_line(r))
    for f in r["failures"]:
        print("  " + f)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(r, f, indent=2)
    sys.exit(1 if r["failures"] else 0)


if __name__ == "__main__":
    main()
