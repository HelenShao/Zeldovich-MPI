#!/usr/bin/env python3
"""Compare the outputs of the seed-noise round-trip test.

Checks (each skipped if its inputs are not given):
  1. noise:  the dumped w = D/sqrt(P) slabs have the expected count/size and <|w|^2> ~ 1
  2. plain:  seed ICs with and without the dump flags are byte identical (dumping does not touch the RNG)
  3. D:      D slabs of the seed run vs the LOAD_EXTERNAL_NOISE rerun
  4. ICs:    ic_* / dens_* files of the seed run vs the rerun (RVZel records, float32 density)

Pass criteria (plan): D and IC floats bitwise identical for > 99.9 % of values, max relative
difference ~1.2e-7 (1 float ulp) for D, and max |diff| / rms at the 1e-6 level for IC fields.
"""
import argparse
import os
import sys

import numpy as np

RVZEL = np.dtype([("ijk", "<u2", (4,)), ("f", "<f4", (6,))])
FIELDS = ("dx", "dy", "dz", "vx", "vy", "vz")


def list_files(root):
    out = []
    for dirpath, _, files in os.walk(root, followlinks=True):
        for f in files:
            out.append(os.path.relpath(os.path.join(dirpath, f), root))
    return sorted(out)


def float_stats(a, b):
    """Bitwise-identical fraction, max |diff|, rms(a), max ulp distance for float arrays."""
    a = np.ascontiguousarray(a).ravel()
    b = np.ascontiguousarray(b).ravel()
    itype = np.int32 if a.dtype == np.float32 else np.int64
    ia, ib = a.view(itype).astype(np.int64), b.view(itype).astype(np.int64)
    same = ia == ib
    d = np.abs(a.astype(np.float64) - b.astype(np.float64))
    same_sign = np.sign(a) == np.sign(b)
    ulp = np.where(same_sign, np.abs(ia - ib), np.iinfo(np.int64).max)
    return dict(n=a.size, same=int(same.sum()), maxdiff=float(d.max(initial=0.0)),
                sumsq=float(np.sum(a.astype(np.float64) ** 2)),
                maxulp=int(ulp.max(initial=0)))


def merge(acc, s):
    if acc is None:
        return dict(s)
    acc["n"] += s["n"]
    acc["same"] += s["same"]
    acc["sumsq"] += s["sumsq"]
    acc["maxdiff"] = max(acc["maxdiff"], s["maxdiff"])
    acc["maxulp"] = max(acc["maxulp"], s["maxulp"])
    return acc


def fmt(s):
    rms = np.sqrt(s["sumsq"] / max(s["n"], 1))
    return (f"n={s['n']} identical={s['same'] / max(s['n'], 1):.7f} max|diff|={s['maxdiff']:.3e} "
            f"rms={rms:.3e} max|diff|/rms={s['maxdiff'] / rms if rms > 0 else 0:.3e} max_ulp={s['maxulp']}")


def check_noise(noise_dir, N, lines):
    files = [f for f in os.listdir(noise_dir) if f.startswith("output_k_space_slab_y.")]
    ys = sorted(int(f.rsplit(".", 1)[1]) for f in files)
    ok = ys == list(range(N // 2 + 1))
    sumsq = 0.0
    nnz = 0
    for y in ys:
        path = os.path.join(noise_dir, f"output_k_space_slab_y.{y}")
        if os.path.getsize(path) != N * N * 16:
            ok = False
            lines.append(f"  bad size: {path}")
            continue
        w = np.fromfile(path, dtype=np.complex128)
        nz = w != 0
        nnz += int(nz.sum())
        sumsq += float(np.sum(np.abs(w[nz]) ** 2))
    mean = sumsq / max(nnz, 1)
    tol = 5.0 / np.sqrt(max(nnz, 1))   # |w|^2 is exponential(1): sd of the mean = 1/sqrt(n)
    ok = ok and abs(mean - 1.0) < tol
    lines.append(f"[noise] {len(ys)} slabs (expect {N // 2 + 1}), nonzero modes {nnz}, "
                 f"<|w|^2> = {mean:.6f} (tolerance {tol:.1e}) -> {'PASS' if ok else 'FAIL'}")
    return ok


def check_bytes_identical(ref, test, label, lines):
    fr, ft = list_files(ref), list_files(test)
    if fr != ft:
        lines.append(f"[{label}] file lists differ: {len(fr)} vs {len(ft)} -> FAIL")
        return False
    ndiff = 0
    for f in fr:
        a = np.fromfile(os.path.join(ref, f), dtype=np.uint8)
        b = np.fromfile(os.path.join(test, f), dtype=np.uint8)
        if a.size != b.size or not np.array_equal(a, b):
            ndiff += 1
    ok = ndiff == 0 and len(fr) > 0
    lines.append(f"[{label}] {len(fr)} files, {ndiff} differ -> {'PASS' if ok else 'FAIL'}")
    return ok


def check_D(ref, test, lines):
    fr, ft = list_files(ref), list_files(test)
    if fr != ft or not fr:
        lines.append(f"[D] file lists differ or empty: {len(fr)} vs {len(ft)} -> FAIL")
        return False
    acc = None
    maxrel = 0.0
    for f in fr:
        a = np.fromfile(os.path.join(ref, f), dtype=np.float64)
        b = np.fromfile(os.path.join(test, f), dtype=np.float64)
        acc = merge(acc, float_stats(a, b))
        nz = a != 0
        if nz.any():
            maxrel = max(maxrel, float(np.max(np.abs(a[nz] - b[nz]) / np.abs(a[nz]))))
        if np.any(b[~nz] != 0):
            maxrel = np.inf
    ok = acc["same"] / acc["n"] > 0.999 and maxrel <= 1.2e-7
    lines.append(f"[D] {len(fr)} slabs: {fmt(acc)} max_rel={maxrel:.3e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def canon_ic_files(root):
    """{canonical name: relative path}. Standalone mode-4 output has ic/zNNN/ic_*, embedded has zNNN/ic_*;
    both are keyed without the leading ic/."""
    out = {}
    for rel in list_files(root):
        b = os.path.basename(rel)
        if not (b.startswith("ic_") or b.startswith("dens_")):
            continue
        parts = rel.split(os.sep)
        if parts[0] == "ic" and len(parts) > 1:
            parts = parts[1:]
        out[os.path.join(*parts)] = rel
    return out


def ic_path_for_dens(canon):
    """dens_<n> -> ic_<n> (flat) or dens/zNNN/dens_<n> -> zNNN/ic_<n>_zNNN (mode 4)."""
    d, b = os.path.split(canon)
    n = b[len("dens_"):]
    parts = d.split(os.sep) if d else []
    if parts and parts[0] == "dens":
        parts = parts[1:]
    if parts and parts[-1].startswith("z"):
        return os.path.join(*parts, f"ic_{n}_{parts[-1]}")
    return os.path.join(*parts, "ic_" + n) if parts else "ic_" + n


def ijk_key(ijk):
    return (ijk[:, 0].astype(np.int64) << 32) | (ijk[:, 1].astype(np.int64) << 16) | ijk[:, 2].astype(np.int64)


def check_ics(ref, test, lines, label="IC", fft_tol=None):
    """Bitwise mode (fft_tol None): > 99.9 % identical floats and max |diff| <= 1e-6 rms.
    FFT-roundoff mode: max |diff| <= fft_tol * rms per field (runs with different FFT plans).
    Particles are matched by (i, j, k) when the record order differs (different rank layouts)."""
    cr, ct = canon_ic_files(ref), canon_ic_files(test)
    fr = sorted(cr)
    if fr != sorted(ct) or not fr:
        lines.append(f"[{label}] file lists differ or empty: {len(cr)} vs {len(ct)} -> FAIL")
        if cr and ct:
            only_r, only_t = sorted(set(cr) - set(ct))[:5], sorted(set(ct) - set(cr))[:5]
            lines.append(f"  only in ref: {only_r}  only in test: {only_t}")
        return False
    ok = True
    hdr_bad = 0
    reordered = 0
    perms = {}
    per = {k: None for k in FIELDS}
    dens = None
    nic = ndens = 0
    for f in fr:
        if not os.path.basename(f).startswith("ic_"):
            continue
        pa, pb = os.path.join(ref, cr[f]), os.path.join(test, ct[f])
        a, b = np.fromfile(pa, RVZEL), np.fromfile(pb, RVZEL)
        if a.size != b.size:
            lines.append(f"  particle count differs: {f} ({a.size} vs {b.size})")
            ok = False
            continue
        nic += 1
        if np.any(a["ijk"] != b["ijk"]):
            ka, kb = ijk_key(a["ijk"]), ijk_key(b["ijk"])
            ia, ib = np.argsort(ka, kind="stable"), np.argsort(kb, kind="stable")
            if not np.array_equal(ka[ia], kb[ib]):
                hdr_bad += int(np.sum(ka[ia] != kb[ib]))
                continue
            reordered += 1
            perms[f] = (ia, ib)
            a, b = a[ia], b[ib]
        for i, k in enumerate(FIELDS):
            per[k] = merge(per[k], float_stats(a["f"][:, i], b["f"][:, i]))
    for f in fr:
        if not os.path.basename(f).startswith("dens_"):
            continue
        a = np.fromfile(os.path.join(ref, cr[f]), np.float32)
        b = np.fromfile(os.path.join(test, ct[f]), np.float32)
        if a.size != b.size:
            lines.append(f"  dens size differs: {f}")
            ok = False
            continue
        ndens += 1
        p = perms.get(ic_path_for_dens(f))
        if p is not None:
            a, b = a[p[0]], b[p[1]]
        dens = merge(dens, float_stats(a, b))
    lines.append(f"[{label}] {nic} ic_ files ({reordered} matched by particle index), {ndens} dens files, "
                 f"particle index mismatches: {hdr_bad}")
    ok = ok and hdr_bad == 0 and nic > 0
    tol = 1e-6 if fft_tol is None else fft_tol
    allf = None
    for name, s in [(k, per[k]) for k in FIELDS] + [("dens", dens)]:
        if s is None:
            continue
        lines.append(f"  {name}: {fmt(s)}")
        rms = np.sqrt(s["sumsq"] / s["n"])
        ok = ok and s["maxdiff"] <= tol * rms
        allf = merge(allf, s)
    if allf is not None:
        frac = allf["same"] / allf["n"]
        lines.append(f"  all floats bitwise identical: {frac:.7f}")
        if fft_tol is None:
            ok = ok and frac > 0.999
    crit = "bitwise" if fft_tol is None else f"FFT roundoff, max|diff| <= {fft_tol:g} rms"
    lines.append(f"[{label}] ({crit}) -> {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--N", type=int, required=True)
    p.add_argument("--noise-dir", default="run_seed/noise_slabs")
    p.add_argument("--ic-plain", default="run_plain/ic")
    p.add_argument("--ic-ref", default="run_seed/ic")
    p.add_argument("--ic-test", default="run_load/ic")
    p.add_argument("--d-ref", default="run_seed/D_slabs")
    p.add_argument("--d-test", default="run_load/D_slabs")
    p.add_argument("--fft-tol", type=float, default=None,
                   help="compare ICs to max|diff| <= FFT_TOL * rms instead of bitwise (different FFT plans)")
    p.add_argument("--ic-extra", nargs=3, action="append", default=[], metavar=("LABEL", "REF", "TEST"),
                   help="additional IC comparison (uses --fft-tol)")
    p.add_argument("--out", default="results.txt")
    a = p.parse_args()

    lines = [f"round-trip comparison, N = {a.N}"]
    results = []
    if a.noise_dir and os.path.isdir(a.noise_dir):
        results.append(check_noise(a.noise_dir, a.N, lines))
    if a.ic_plain and os.path.isdir(a.ic_plain):
        results.append(check_bytes_identical(a.ic_plain, a.ic_ref, "plain vs seed+dump ICs", lines))
    if a.d_ref and os.path.isdir(a.d_ref) and a.d_test and os.path.isdir(a.d_test):
        results.append(check_D(a.d_ref, a.d_test, lines))
    results.append(check_ics(a.ic_ref, a.ic_test, lines, "IC", a.fft_tol))
    for label, r, t in a.ic_extra:
        results.append(check_ics(r, t, lines, label, a.fft_tol))
    lines.append(f"OVERALL: {'PASS' if all(results) else 'FAIL'}")
    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
