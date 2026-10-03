#!/usr/bin/env python
"""Particle-by-particle comparison of an Abacus FLAMINGO-matched subsample with a FLAMINGO snapshot.

Particles are matched by Lagrangian lattice index:
  Abacus  : lagr_idx (i, j, k) from the subsample PIDs, axes (x, y, z)
  FLAMINGO: monofonIC DM ParticleIDs = overload * ((i_m * N + j_m) * N + k_m), overload = 2 with baryons
  Abacus (x, y, z) = monofonIC (z, y, x), so (i, j, k)_Abacus = (k_m, j_m, i_m).

Reports (to --out, a text file): match count, periodic position offsets, displacement and velocity
correlations (sign check: expect r ~ +1), and the density cross-correlation r(k) and P ratio of the
matched set. Optionally a FRESCO-style figure (--plot).

Run on a compute node (ASDF and multi-GB HDF5 reads), e.g. via pbs_compare.sh.
"""
import argparse
import glob
import os
import sys
from pathlib import Path

import h5py
import numpy as np


class Tee:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(s + "\n")
        self.f.flush()


def load_abacus(zdir, log):
    from abacusnbody.data.read_abacus import read_asdf

    files = sorted(glob.glob(os.path.join(zdir, "col*", "particles_*.asdf")))
    if not files:
        sys.exit(f"no subsample files under {zdir}")
    pos, vel, lagr = [], [], []
    header = None
    for fn in files:
        # output-particle format: pid = (N, 3) uint16 Lagrangian lattice index, undefined for MAP rows
        t = read_asdf(fn, load=("pos", "vel", "pid", "is_map"), verbose=False)
        header = header or dict(t.meta)
        dm = ~np.asarray(t["is_map"], dtype=bool)
        pos.append(np.asarray(t["pos"], dtype=np.float64)[dm])
        vel.append(np.asarray(t["vel"], dtype=np.float32)[dm])
        lagr.append(np.asarray(t["pid"], dtype=np.int64)[dm])
    pos, vel, lagr = np.concatenate(pos), np.concatenate(vel), np.concatenate(lagr)
    L, N = float(header["BoxSize"]), int(round(header["NP"] ** (1 / 3)))
    log(f"[abacus] {len(files)} files, {len(pos):,} particles ({len(pos) / N**3:.4%} of {N}^3), "
        f"L = {L} Mpc/h, z = {header.get('Redshift')}")
    return pos, vel, lagr, L, N, header


def flamingo_files(path):
    p = Path(path)
    files = sorted(p.glob("*.hdf5")) if p.is_dir() else [p]
    if not files:
        sys.exit(f"no FLAMINGO HDF5 files at {path}")
    return files


def phys_factor(dset):
    """Conversion of a SWIFT dataset to physical cgs, cosmological factors included."""
    for key in ("Conversion factor to physical CGS (including cosmological corrections)",):
        if key in dset.attrs:
            return float(np.squeeze(dset.attrs[key]))
    raise KeyError(f"{dset.name}: no physical CGS conversion attribute")


def load_flamingo_matched(files, target_ids, chunk, log):
    """Stream PartType1 and keep particles whose ID is in target_ids (sorted int64)."""
    with h5py.File(files[0], "r") as f:
        hdr = dict(f["Header"].attrs)
        cos = dict(f["Cosmology"].attrs) if "Cosmology" in f else {}
    z = float(np.squeeze(hdr.get("Redshift")))
    a = 1.0 / (1.0 + z)
    h = float(np.squeeze(cos.get("h", np.nan)))
    out_id, out_pos, out_vel = [], [], []
    ntot = 0
    for fn in files:
        with h5py.File(fn, "r") as f:
            ids_d, pos_d, vel_d = f["PartType1/ParticleIDs"], f["PartType1/Coordinates"], f["PartType1/Velocities"]
            cm_per_mpc = 3.0856775814913673e24
            pos_fac = phys_factor(pos_d) / cm_per_mpc / a * h  # -> comoving Mpc/h
            vel_fac = phys_factor(vel_d) / 1e5                 # -> peculiar km/s
            n = ids_d.shape[0]
            for s in range(0, n, chunk):
                ids = ids_d[s:s + chunk].astype(np.int64)
                j = np.searchsorted(target_ids, ids)
                j[j == len(target_ids)] = 0
                keep = target_ids[j] == ids
                if keep.any():
                    out_id.append(ids[keep])
                    out_pos.append(pos_d[s:s + chunk][keep].astype(np.float64) * pos_fac)
                    out_vel.append(vel_d[s:s + chunk][keep].astype(np.float32) * vel_fac)
            ntot += n
    log(f"[flamingo] {len(files)} file(s), {ntot:,} DM particles scanned, z = {z:.4f}, h = {h}")
    return np.concatenate(out_id), np.concatenate(out_pos), np.concatenate(out_vel), hdr, z, h


def periodic(d, L):
    return d - L * np.round(d / L)


def cic_k(pos, L, ng):
    """CIC density contrast in k-space (rfftn), positions in [0, L)."""
    x = (pos / L * ng) % ng
    i0 = np.floor(x).astype(np.int64)
    w1 = x - i0
    rho = np.zeros(ng**3)
    for dx in (0, 1):
        wx = w1[:, 0] if dx else 1 - w1[:, 0]
        for dy in (0, 1):
            wy = w1[:, 1] if dy else 1 - w1[:, 1]
            for dz in (0, 1):
                wz = w1[:, 2] if dz else 1 - w1[:, 2]
                flat = (((i0[:, 0] + dx) % ng) * ng + (i0[:, 1] + dy) % ng) * ng + (i0[:, 2] + dz) % ng
                rho += np.bincount(flat, weights=wx * wy * wz, minlength=ng**3)
    rho = rho.reshape(ng, ng, ng)
    rho = rho / rho.mean() - 1
    return np.fft.rfftn(rho)


def kbins(ng, L, nb):
    kf = 2 * np.pi / L
    kx = np.fft.fftfreq(ng, 1 / ng)
    kz = np.fft.rfftfreq(ng, 1 / ng)
    kk = np.sqrt(kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kz[None, None, :] ** 2)
    edges = np.logspace(0, np.log10(ng / 2), nb + 1)
    idx = np.digitize(kk, edges) - 1
    return idx, edges * kf


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abacus", required=True, help="subsample redshift dir, e.g. .../subsamples/z0.500")
    ap.add_argument("--flamingo", help="FLAMINGO snapshot: a .hdf5 (virtual or single) or a dir of chunks")
    ap.add_argument("--abacus-only", action="store_true", help="load and summarise the Abacus subsample, then stop")
    ap.add_argument("--overload", type=int, default=0, help="ID multiplier: 1 for DMO, 2 with baryons (0 = infer)")
    ap.add_argument("--chunk", type=int, default=50_000_000)
    ap.add_argument("--ngrid", type=int, default=256, help="mesh for r(k) of the matched set")
    ap.add_argument("--out", default="compare_results.txt")
    ap.add_argument("--plot", help="optional figure path (.png/.pdf)")
    args = ap.parse_args()
    log = Tee(args.out)
    log(f"abacus: {args.abacus}\nflamingo: {args.flamingo}")

    pos_a, vel_a, lagr, L_a, N, hdr_a = load_abacus(args.abacus, log)
    # Abacus (i, j, k) -> monofonIC (i_m, j_m, k_m) = (k, j, i)
    cell = (lagr[:, 2] * N + lagr[:, 1]) * N + lagr[:, 0]
    if args.abacus_only:
        log(f"[abacus] header keys: {sorted(hdr_a)[:40]}")
        log(f"[abacus] pos range {pos_a.min(axis=0)} .. {pos_a.max(axis=0)}; lagr_idx range "
            f"{lagr.min(axis=0)} .. {lagr.max(axis=0)}; unique cells {len(np.unique(cell)):,}")
        log(f"[abacus] vel rms per axis [km/s] {vel_a.std(axis=0)}")
        return
    if not args.flamingo:
        sys.exit("--flamingo is required unless --abacus-only")

    files = flamingo_files(args.flamingo)
    overload = args.overload
    if overload == 0:
        with h5py.File(files[0], "r") as f:
            ngas = int(np.atleast_1d(f["Header"].attrs.get("NumPart_Total", [0]))[0])
        overload = 2 if ngas > 0 else 1
    log(f"ID overload = {overload} (DM ID = overload * cell)")
    target = overload * cell
    order = np.argsort(target)
    target_sorted = target[order]

    ids_f, pos_f, vel_f, hdr_f, z_f, h_f = load_flamingo_matched(files, target_sorted, args.chunk, log)
    L_f = float(np.atleast_1d(hdr_f["BoxSize"])[0]) * h_f  # comoving Mpc -> Mpc/h
    log(f"[flamingo] L = {L_f:.4f} Mpc/h; matched {len(ids_f):,} of {len(target):,} Abacus particles "
        f"({len(ids_f) / len(target):.4%})")
    if len(ids_f) == 0:
        sys.exit("no matches: check --overload and the axis/ID convention")

    # align: FLAMINGO order -> Abacus rows
    rows = order[np.searchsorted(target_sorted, ids_f)]
    pa, va = pos_a[rows], vel_a[rows]
    # FLAMINGO (x, y, z) -> Abacus axes: swap x and z; compare in box fractions scaled to L_a
    pf = pos_f[:, ::-1] / L_f * L_a
    vf = vel_f[:, ::-1]
    pa0 = pa + L_a / 2  # Abacus box is centred on 0

    d = periodic(pa0 - pf, L_a)
    off = np.median(d, axis=0)
    log(f"median offset Abacus - FLAMINGO [Mpc/h] = {off}")
    dr = np.linalg.norm(periodic(d - off, L_a), axis=1)
    pct = np.percentile(dr, [16, 50, 84, 95, 99])
    log(f"|dx| after removing the median offset [Mpc/h]: p16 {pct[0]:.3f}  p50 {pct[1]:.3f}  p84 {pct[2]:.3f}  "
        f"p95 {pct[3]:.3f}  p99 {pct[4]:.3f}")
    for thr in (0.1, 0.25, 0.5, 1.0, 2.0):
        log(f"  fraction with |dx| < {thr:4.2f} Mpc/h: {np.mean(dr < thr):.4f}")

    # displacements from the shared lattice (monofonIC q = (i_m, j_m, k_m) / N, in Abacus axes = lagr / N)
    q = lagr[rows] / N * L_a
    sa = periodic(pa0 - q, L_a)
    sa -= np.median(sa, axis=0)
    sf = periodic(pf + off - q, L_a)
    sf -= np.median(sf, axis=0)
    for c, name in enumerate("xyz"):
        r_s = np.corrcoef(sa[:, c], sf[:, c])[0, 1]
        r_v = np.corrcoef(va[:, c], vf[:, c])[0, 1]
        log(f"axis {name}: corr(displacement) = {r_s:+.5f}   corr(velocity) = {r_v:+.5f}   "
            f"rms disp A/F = {sa[:, c].std():.3f}/{sf[:, c].std():.3f}   rms vel A/F = {va[:, c].std():.1f}/{vf[:, c].std():.1f}")
    log("(sign check: correlations ~ +1 mean same-sign ICs; ~ -1 would mean inverted phases)")

    # density of the matched set: cross-correlation and power ratio
    ng = args.ngrid
    da = cic_k(pa0 % L_a, L_a, ng)
    df = cic_k((pf + off) % L_a, L_a, ng)
    idx, kedges = kbins(ng, L_a, 16)
    log(f"density of the matched set on {ng}^3 (k in h/Mpc):   k_lo     k_hi     r(k)   P_A/P_F   N_modes")
    rows_out = []
    for b in range(len(kedges) - 1):
        m = idx == b
        if not m.any():
            continue
        paa = np.mean(np.abs(da[m]) ** 2)
        pff = np.mean(np.abs(df[m]) ** 2)
        pxf = np.mean((da[m] * np.conj(df[m])).real)
        r = pxf / np.sqrt(paa * pff)
        rows_out.append((kedges[b], kedges[b + 1], r, paa / pff))
        log(f"   {kedges[b]:8.4f} {kedges[b + 1]:8.4f}  {r:+.5f}  {paa / pff:8.5f}  {m.sum():8d}")

    if args.plot:
        sys.path.insert(0, "/home/helenshao/InitialConditions/hackAurora")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from plot_style import apply_plot_style, palette_for

        apply_plot_style("fresco")
        pal = palette_for("fresco")
        arr = np.array(rows_out)
        kc = np.sqrt(arr[:, 0] * arr[:, 1])
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        ax[0].semilogx(kc, arr[:, 2], "o-", color=pal["blue"], label="r(k)")
        ax[0].semilogx(kc, arr[:, 3], "s--", color=pal["orange"], label=r"$P_{\rm Abacus}/P_{\rm FLAMINGO}$")
        ax[0].axhline(1, color="k", lw=0.6)
        ax[0].set_xlabel(r"$k\ [h\,{\rm Mpc}^{-1}]$")
        ax[0].legend()
        ax[1].hist(dr, bins=np.logspace(-3, 1.5, 80), color=pal["blue"], histtype="step")
        ax[1].set_xscale("log")
        ax[1].set_xlabel(r"$|\Delta x|\ [h^{-1}{\rm Mpc}]$")
        ax[1].set_ylabel("particles")
        fig.suptitle(f"cosm202 vs FLAMINGO, z = {z_f:.2f}")
        fig.tight_layout()
        fig.savefig(args.plot)
        log(f"figure: {os.path.abspath(args.plot)}")

    log(f"results written to {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
