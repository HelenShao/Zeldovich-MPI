#!/usr/bin/env python3
"""Regression tests: C Kernel_MPI against the confirmed kernel_MPI.py.

Tiers
  quick  (default, seconds to minutes) synthetic inputs, N in {16, 32}, all four f_NL types,
         several rank layouts. Per case:
           [py]     C (compat, M=2N) vs kernel_MPI.py, same ranks          -> must pass
           [direct] C vs reference_direct.py (plain numpy, independent)     -> must pass
         Per (N, f_NL) group:
           [ranks]  C at every rank layout vs C on 1 rank                   -> must pass
                    forced staged (int-overflow-safe) Alltoallv vs normal   -> must pass
           [nyq]    C with KERNEL_ZERO_NYQUIST=1 vs direct numpy, same zeroing -> must pass
                    ... and vs compat mode (size of the Nyquist effect)     -> report only
           [pad]    zero-Nyquist C with M = 3N/2 vs M = 2N (alias-free)     -> must pass
                    compat C with M = 3N/2 vs M = 2N                        -> report only
  full   the real dumps in 1rank-run-cosm500/bin_files (or --full-bin-dir), all four
         f_NL types, C vs kernel_MPI.py with the dump's rank count. Needs a machine
         that can hold kernel_MPI.py's (2N)^3 arrays (~100+ GB at N=500).

Python reference outputs are cached in the work directory, keyed by a hash of the
inputs and of kernel_MPI.py, so they are only recomputed when either changes.

Exit status: 0 all pass, 1 any failure, 77 skipped (missing pyfftw / mpiexec / binary).
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

from compare_kernel_outputs import compare, summary_line  # noqa: E402
from make_test_phi import write_inputs, write_par2  # noqa: E402

KERNEL_PY = os.path.join(REPO, "1rank-run-cosm500", "kernel_MPI.py")
FNL_TYPES = ["local", "equil", "orthoCMB", "orthoLSS"]
BOX_SIZE = 500.0
SEED = 20261017

QUICK_RANK_LAYOUTS = [(1, 1), (2, 1), (2, 2), (3, 1), (4, 2), (6, 2)]   # (ranks, NumZRanks)
QUICK_N = {16: QUICK_RANK_LAYOUTS, 32: [(1, 1), (3, 1), (4, 2)]}


def log(msg):
    print(msg, flush=True)


def choose_3half_M(N):
    M = (3 * N + 1) // 2
    while M % 4 or (M - N) % 2:
        M += 1
    return M


class Runner:
    def __init__(self, args):
        self.args = args
        self.work = os.path.abspath(args.work_dir)
        self.exe = os.path.abspath(args.kernel_exe or os.path.join(args.build_dir, "src", "Kernel_MPI"))
        self.results = []
        self.mpiexec = args.mpiexec.split()
        if args.oversubscribe:
            self.mpiexec.append("--oversubscribe")

    # ---------------------------------------------------------------- running things
    def run(self, cmd, env_extra, log_path, cwd=None):
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = str(self.args.omp_threads)
        env.update(env_extra)
        t0 = time.perf_counter()
        with open(log_path, "w") as f:
            p = subprocess.run(cmd, env=env, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
        if p.returncode != 0:
            tail = open(log_path).read()[-3000:]
            raise RuntimeError(f"command failed ({p.returncode}): {' '.join(cmd)}\n--- {log_path} ---\n{tail}")
        return time.perf_counter() - t0

    def run_c(self, bin_dir, par, ranks, out_dir, env_extra=None):
        shutil.rmtree(out_dir, ignore_errors=True)
        os.makedirs(out_dir)
        cmd = self.mpiexec + ["-n", str(ranks), self.exe, par, "--bin-dir", bin_dir, "--out-dir", out_dir]
        return self.run(cmd, env_extra or {}, out_dir + ".log")

    def run_python_ref(self, bin_dir, par, ranks, out_dir):
        """kernel_MPI.py writes to $BIN_DIR/bin_files_output; run it on a private copy."""
        key = hashlib.sha256()
        for path in sorted(os.listdir(bin_dir)):
            if path.endswith(".bin"):
                with open(os.path.join(bin_dir, path), "rb") as f:
                    key.update(path.encode())
                    key.update(hashlib.sha256(f.read()).digest())
        key.update(open(par, "rb").read())
        key.update(open(KERNEL_PY, "rb").read())
        key.update(str(ranks).encode())
        stamp = os.path.join(out_dir, ".hash")
        if os.path.exists(stamp) and open(stamp).read() == key.hexdigest():
            return 0.0

        run_dir = out_dir + "_run"
        shutil.rmtree(run_dir, ignore_errors=True)
        os.makedirs(run_dir)
        for name in os.listdir(bin_dir):
            if name.startswith("Phi_gaussian_rank") and name.endswith(".bin"):
                os.symlink(os.path.join(os.path.abspath(bin_dir), name), os.path.join(run_dir, name))
        env = {"PARAM_FILE": os.path.abspath(par), "BIN_DIR": run_dir, "NUM_RANKS": str(ranks)}
        dt = self.run([self.args.python, KERNEL_PY], env, out_dir + ".log", cwd=os.path.dirname(KERNEL_PY))
        shutil.rmtree(out_dir, ignore_errors=True)
        shutil.move(os.path.join(run_dir, "bin_files_output"), out_dir)
        shutil.rmtree(run_dir)
        # kernel_MPI.py names its planes phi_py_output_*; the comparer expects the C name.
        for name in os.listdir(out_dir):
            if name.startswith("phi_py_output_"):
                os.rename(os.path.join(out_dir, name),
                          os.path.join(out_dir, name.replace("phi_py_output_", "phi_output_", 1)))
        with open(stamp, "w") as f:
            f.write(key.hexdigest())
        return dt

    def run_direct_ref(self, bin_dir, ranks, fnl, out_dir, zero_nyquist=False):
        shutil.rmtree(out_dir, ignore_errors=True)
        cmd = [self.args.python, os.path.join(HERE, "reference_direct.py"), "--bin-dir", bin_dir,
               "--out-dir", out_dir, "--ranks", str(ranks), "--box-size", str(BOX_SIZE), "--fnl", fnl]
        if zero_nyquist:
            cmd.append("--zero-nyquist")
        return self.run(cmd, {}, out_dir + ".log")

    # ---------------------------------------------------------------- bookkeeping
    def check(self, label, kind, dir_a, dir_b, gate=True, rel_l2_tol=1e-10, max_rel_tol=1e-8):
        r = compare(dir_a, dir_b, rel_l2_tol, max_rel_tol)
        ok = not r["failures"] or not gate
        status = ("PASS" if not r["failures"] else "FAIL") if gate else "INFO"
        log(f"  {status:4s} [{kind:6s}] {label:34s} {summary_line(r)}")
        if r["failures"] and gate:
            for f in r["failures"]:
                log(f"         {f}")
        r.update({"label": label, "kind": kind, "status": status})
        self.results.append(r)
        return ok

    # ---------------------------------------------------------------- tiers
    def quick(self):
        root = os.path.join(self.work, "quick")
        for N, layouts in QUICK_N.items():
            if self.args.N and N not in self.args.N:
                continue
            for fnl in FNL_TYPES:
                if self.args.fnl and fnl not in self.args.fnl:
                    continue
                log(f"\n=== N={N} f_NL={fnl} ===")
                c_dirs = {}
                for ranks, nz in layouts:
                    tag = f"N{N}_R{ranks}_nz{nz}"
                    case = os.path.join(root, f"N{N}", fnl, tag)
                    bin_dir = os.path.join(root, f"N{N}", f"inputs_R{ranks}")
                    if not os.path.isdir(bin_dir):
                        write_inputs(bin_dir, N, ranks, SEED)
                    os.makedirs(case, exist_ok=True)
                    par = os.path.join(case, "param.par2")
                    write_par2(par, BOX_SIZE, nz, fnl)

                    c_out = os.path.join(case, "c_out")
                    self.run_c(bin_dir, par, ranks, c_out)
                    c_dirs[(ranks, nz)] = (c_out, bin_dir, par)

                    py_out = os.path.join(case, "py_out")
                    self.run_python_ref(bin_dir, par, ranks, py_out)
                    self.check(f"{tag} C vs kernel_MPI.py", "py", c_out, py_out)

                    direct_out = os.path.join(case, "direct_out")
                    self.run_direct_ref(bin_dir, ranks, fnl, direct_out)
                    self.check(f"{tag} C vs direct numpy", "direct", c_out, direct_out)

                # Rank invariance: every layout vs the first (1 rank)
                base = layouts[0]
                for lay in layouts[1:]:
                    self.check(f"N{N} R{lay[0]}nz{lay[1]} vs R{base[0]}nz{base[1]}", "ranks",
                               c_dirs[lay][0], c_dirs[base][0])

                # Overflow-safe (staged) Alltoallv path, normally only used past INT_MAX offsets
                lay = layouts[-1]
                staged_out = c_dirs[lay][0] + "_staged_a2a"
                self.run_c(c_dirs[lay][1], c_dirs[lay][2], lay[0], staged_out, {"KERNEL_FORCE_STAGED_A2A": "1"})
                self.check(f"N{N} R{lay[0]} staged Alltoallv", "ranks", staged_out, c_dirs[lay][0])

                # Nyquist zeroing: C vs the direct reference with the same zeroing (gated),
                # and the size of the change relative to compat mode (informational)
                lay = layouts[-1]
                c_out, bin_dir, par = c_dirs[lay]
                nyq_out = c_out + "_zero_nyquist"
                self.run_c(bin_dir, par, lay[0], nyq_out, {"KERNEL_ZERO_NYQUIST": "1"})
                direct_nyq = os.path.join(os.path.dirname(c_out), "direct_out_zero_nyquist")
                self.run_direct_ref(bin_dir, lay[0], fnl, direct_nyq, zero_nyquist=True)
                self.check(f"N{N} zero-Nyquist C vs direct numpy", "nyq", nyq_out, direct_nyq)
                self.check(f"N{N} zero-Nyquist vs compat (size)", "nyq", nyq_out, c_out, gate=False)

                # Padding: M = 3N/2 vs M = 2N. Exact (round-off) once the Nyquist modes are
                # zeroed; in compat mode the unpaired -N/2 modes make the packed product
                # M-dependent, so that difference is reported but not gated.
                M32 = choose_3half_M(N)
                pad_nyq = nyq_out + f"_M{M32}"
                self.run_c(bin_dir, par, lay[0], pad_nyq, {"KERNEL_ZERO_NYQUIST": "1", "KERNEL_PAD_M": str(M32)})
                self.check(f"N{N} zero-Nyq M={M32} vs M={2 * N}", "pad", pad_nyq, nyq_out)
                pad_out = c_out + f"_M{M32}"
                self.run_c(bin_dir, par, lay[0], pad_out, {"KERNEL_PAD_M": str(M32)})
                self.check(f"N{N} compat M={M32} vs M={2 * N}", "pad", pad_out, c_out, gate=False)

    def full(self):
        bin_dir = os.path.abspath(self.args.full_bin_dir)
        names = [n for n in os.listdir(bin_dir) if n.startswith("Phi_gaussian_rank") and n.endswith(".bin")]
        if not names:
            raise SystemExit(f"no Phi_gaussian_rank*_y*.bin in {bin_dir}")
        ranks = max(int(n.split("rank")[1].split("_")[0]) for n in names) + 1
        nz = self.args.full_num_z_ranks
        root = os.path.join(self.work, "full")
        for fnl in FNL_TYPES:
            if self.args.fnl and fnl not in self.args.fnl:
                continue
            log(f"\n=== full: {bin_dir} R={ranks} NumZRanks={nz} f_NL={fnl} ===")
            case = os.path.join(root, fnl)
            os.makedirs(case, exist_ok=True)
            par = os.path.join(case, "param.par2")
            write_par2(par, self.args.full_box_size, nz, fnl)
            c_out = os.path.join(case, "c_out")
            dt_c = self.run_c(bin_dir, par, ranks, c_out)
            py_out = os.path.join(case, "py_out")
            dt_py = self.run_python_ref(bin_dir, par, ranks, py_out)
            log(f"  C {dt_c:.1f} s, kernel_MPI.py {dt_py:.1f} s (0 = cached)")
            self.check(f"full {fnl} C vs kernel_MPI.py", "py", c_out, py_out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["quick", "full"], default="quick")
    ap.add_argument("--build-dir", default=os.path.join(REPO, "build", "kernel_double"))
    ap.add_argument("--kernel-exe", help="path to Kernel_MPI (default <build-dir>/src/Kernel_MPI)")
    ap.add_argument("--work-dir", default=os.path.join(REPO, "build", "kernel_regression_work"))
    ap.add_argument("--python", default=sys.executable, help="interpreter with numpy + pyfftw")
    ap.add_argument("--mpiexec", default="mpiexec")
    ap.add_argument("--oversubscribe", action="store_true", help="pass --oversubscribe (Open MPI)")
    ap.add_argument("--omp-threads", type=int, default=1)
    ap.add_argument("--N", type=int, action="append", help="quick tier: restrict to these N")
    ap.add_argument("--fnl", action="append", choices=FNL_TYPES, help="restrict to these f_NL types")
    ap.add_argument("--full-bin-dir", default=os.path.join(REPO, "1rank-run-cosm500", "bin_files"))
    ap.add_argument("--full-box-size", type=float, default=None,
                    help="full tier BoxSize (default: read from 1rank-run-cosm500/.runtime/param.par2)")
    ap.add_argument("--full-num-z-ranks", type=int, default=None)
    ap.add_argument("--json", help="write all comparison results here")
    args = ap.parse_args()

    runner = Runner(args)
    if not os.access(runner.exe, os.X_OK):
        log(f"SKIP: {runner.exe} not built (run tests/kernel_regression/build_kernel.sh)")
        sys.exit(77)
    if shutil.which(runner.mpiexec[0]) is None:
        log(f"SKIP: {runner.mpiexec[0]} not found")
        sys.exit(77)
    if subprocess.run([args.python, "-c", "import pyfftw, numpy"], capture_output=True).returncode != 0:
        log(f"SKIP: {args.python} cannot import pyfftw/numpy (needed to run kernel_MPI.py)")
        sys.exit(77)

    if args.tier == "full":
        runtime_par = os.path.join(REPO, "1rank-run-cosm500", ".runtime", "param.par2")
        if args.full_box_size is None or args.full_num_z_ranks is None:
            sys.path.insert(0, os.path.dirname(KERNEL_PY))
            vals = {}
            if os.path.exists(runtime_par):
                for line in open(runtime_par):
                    line = line.split("#", 1)[0]
                    if "=" in line:
                        k, v = line.split("=", 1)
                        vals[k.strip()] = v.strip().strip('"').strip("'")
            if args.full_box_size is None:
                if "BoxSize" not in vals:
                    raise SystemExit("--full-box-size not given and no BoxSize in " + runtime_par)
                args.full_box_size = float(vals["BoxSize"])
            if args.full_num_z_ranks is None:
                args.full_num_z_ranks = int(float(vals.get("NumZRanks", 1)))

    os.makedirs(args.work_dir, exist_ok=True)
    log(f"Kernel_MPI regression ({args.tier}): exe={runner.exe}\n  reference={KERNEL_PY}\n  work={runner.work}")
    t0 = time.perf_counter()
    try:
        runner.quick() if args.tier == "quick" else runner.full()
    except RuntimeError as e:
        log(f"\nERROR: {e}")
        sys.exit(1)

    gated = [r for r in runner.results if r["status"] != "INFO"]
    failed = [r for r in gated if r["status"] == "FAIL"]
    log(f"\n{len(gated) - len(failed)}/{len(gated)} checks passed "
        f"({sum(r['status'] == 'INFO' for r in runner.results)} informational) in {time.perf_counter() - t0:.1f} s")
    worst = max((r.get("rel_l2", np.inf) for r in gated), default=0.0)
    log(f"worst gated rel_l2 = {worst:.3e}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(runner.results, f, indent=2)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
