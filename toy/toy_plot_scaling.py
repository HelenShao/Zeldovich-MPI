#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

# Threads
threads = np.array([1, 2, 4, 8, 16], dtype=int)

# Case A: No MPI (single process)
zloop_no_mpi = np.array([1.012131, 0.547590, 0.305208, 0.162468, 0.091175])

# Case B: 9 MPI ranks (approx mean per-rank times)
zloop_mpi9 = np.array([0.803, 0.816, 1.034, 1.067, 1.361])

# Case C: 1 MPI rank (mpiexec -n 1)
zloop_mpi1 = np.array([0.804461, 0.813792, 1.002895, 1.069554, 1.257730])

# Case D: 1 MPI rank + TOY_FIRST_TOUCH=1
zloop_mpi1_ft = np.array([0.803098, 0.819353, 0.993355, 1.098941, 1.322112])

cases = {
    "No MPI (single process)": zloop_no_mpi,
    "MPI 9 ranks (per rank)": zloop_mpi9,
    "MPI 1 rank": zloop_mpi1,
    "MPI 1 rank + first-touch": zloop_mpi1_ft,
}

colors = {
    "No MPI (single process)": "C0",
    "MPI 9 ranks (per rank)": "C1",
    "MPI 1 rank": "C2",
    "MPI 1 rank + first-touch": "C3",
}

# --- 1) Line plot: Z-loop time vs threads -----------------------------------
plt.figure(figsize=(6, 4))
for name, vals in cases.items():
    plt.plot(threads, vals, marker="o", label=name, color=colors[name])
plt.xticks(threads)
plt.xlabel("OMP threads")
plt.ylabel("Z-loop time [s] (10 Y-repeats)")
plt.title("Toy z-loop: Z-loop time vs OMP threads")
plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("toy_zloop_time_vs_threads.png", dpi=200)

# --- 2) Line plot: Speedup vs threads (per case) ----------------------------
plt.figure(figsize=(6, 4))
for name, vals in cases.items():
    speedup = vals[0] / vals  # relative to 1 thread in that case
    plt.plot(threads, speedup, marker="o", label=name, color=colors[name])
plt.xticks(threads)
plt.xlabel("OMP threads")
plt.ylabel("Speedup (T1 / Tn)")
plt.title("Toy z-loop: Speedup vs OMP threads")
plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("toy_zloop_speedup_vs_threads.png", dpi=200)

# --- 3) Histograms of Z-loop times (simple view) ----------------------------
# Flatten times for each case into separate subplot
fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)
axes = axes.ravel()
bins = 8

for ax, (name, vals) in zip(axes, cases.items()):
    ax.bar(threads, vals, width=0.6, color=colors[name])
    ax.set_title(name, fontsize=9)
    ax.set_xlabel("OMP threads")
    ax.set_ylabel("Z-loop time [s]")
    ax.set_xticks(threads)
    ax.grid(True, axis="y", linestyle="--", alpha=0.3)

fig.suptitle("Toy z-loop: Z-loop time per case", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("toy_zloop_time_histograms.png", dpi=200)

print("Wrote plots:")
print("  toy_zloop_time_vs_threads.png")
print("  toy_zloop_speedup_vs_threads.png")
print("  toy_zloop_time_histograms.png")