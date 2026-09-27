#!/usr/bin/env python3
"""
Parse OMP sweep timing summaries under runs/omp-*/ and generate scaling plots.

Ignores legacy "Stage 2 (Metadata exchange)" when unused.
Maps Stage 3 (Communication) -> Stage 2, Stage 4 (Streaming) -> Stage 3.

Outputs (under --results):
- omp_stage_times.png
- omp_stage_histogram.png
- omp_stage_speedup.png
- omp_stage_efficiency.png
- omp_stage_summary.csv

Usage:
  python3 plot_omp_stage_scaling.py
  python3 plot_omp_stage_scaling.py --results results
  python3 plot_omp_stage_scaling.py --base-dir /path/to/experiment
  python3 plot_omp_stage_scaling.py --runs-dir /path/to/experiment/runs --results /path/to/experiment/results
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


STAGE1_RE = re.compile(r"Stage 1 \(Y-slice generation \+ 2D FFT\):\s*([\d.]+)\s*s")
COMM_RE = re.compile(r"Stage 3 \(Communication: Alltoallv\):\s*([\d.]+)\s*s")
STREAM_RE = re.compile(r"Stage 4 \(Streaming: Unpack\+FFT\+Write\):\s*([\d.]+)\s*s")
TOTAL_RE = re.compile(r"Total time \(including all stages\):\s*([\d.]+)\s*s")


def parse_log(path: Path) -> dict[str, float]:
    text = path.read_text()
    stage1 = STAGE1_RE.search(text)
    stage2 = COMM_RE.search(text)
    stage3 = STREAM_RE.search(text)
    total = TOTAL_RE.search(text)

    if not (stage1 and stage2 and stage3 and total):
        raise ValueError(f"Could not parse timing summary from {path}")

    return {
        "Stage 1": float(stage1.group(1)),
        "Stage 2": float(stage2.group(1)),
        "Stage 3": float(stage3.group(1)),
        "Total": float(total.group(1)),
    }


def collect_rows(logs: list[tuple[int, Path]]) -> list[tuple[int, dict[str, float]]]:
    rows: list[tuple[int, dict[str, float]]] = []
    for threads, path in logs:
        try:
            rows.append((threads, parse_log(path)))
        except ValueError:
            print(f"Skipping incomplete or unparsable log: {path}", file=sys.stderr)
    return rows


def resolve_default_base_dir() -> Path:
    script_dir = Path(__file__).resolve().parent
    experiment_root = script_dir.parent
    if (experiment_root / "runs").is_dir():
        return experiment_root
    return script_dir


def resolve_run_log(omp_dir: Path, threads: int) -> Path | None:
    matches = sorted(omp_dir.glob(f"*OMP{threads}_bin_generation.log"))
    if matches:
        return matches[0]

    candidates = sorted(omp_dir.glob("*_bin_generation.log"))
    if not candidates:
        return None

    needle = f"OMP{threads}_"
    for candidate in candidates:
        if needle in candidate.name:
            return candidate

    if len(candidates) == 1:
        return candidates[0]

    return None


def discover_logs(runs_dir: Path) -> list[tuple[int, Path]]:
    logs: list[tuple[int, Path]] = []

    for omp_dir in sorted(runs_dir.glob("omp-*")):
        if not omp_dir.is_dir():
            continue

        match = re.fullmatch(r"omp-(\d+)", omp_dir.name)
        if not match:
            continue

        thread_count = int(match.group(1))
        run_log = resolve_run_log(omp_dir, thread_count)
        if run_log is not None:
            logs.append((thread_count, run_log))

    logs.sort(key=lambda item: item[0])
    return logs


def write_csv(rows: list[tuple[int, dict[str, float]]], out_path: Path) -> None:
    baseline_threads = rows[0][0]
    baseline = rows[0][1]
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "omp_threads",
                "baseline_omp_threads",
                "stage_1_time_s",
                "stage_2_time_s",
                "stage_3_time_s",
                "total_time_s",
                "stage_1_speedup",
                "stage_2_speedup",
                "stage_3_speedup",
                "total_speedup",
                "stage_1_efficiency",
                "stage_2_efficiency",
                "stage_3_efficiency",
                "total_efficiency",
            ]
        )

        for threads, data in rows:
            thread_scale = threads / baseline_threads
            writer.writerow(
                [
                    threads,
                    baseline_threads,
                    data["Stage 1"],
                    data["Stage 2"],
                    data["Stage 3"],
                    data["Total"],
                    baseline["Stage 1"] / data["Stage 1"],
                    baseline["Stage 2"] / data["Stage 2"],
                    baseline["Stage 3"] / data["Stage 3"],
                    baseline["Total"] / data["Total"],
                    ((baseline["Stage 1"] / data["Stage 1"]) / thread_scale) * 100.0,
                    ((baseline["Stage 2"] / data["Stage 2"]) / thread_scale) * 100.0,
                    ((baseline["Stage 3"] / data["Stage 3"]) / thread_scale) * 100.0,
                    ((baseline["Total"] / data["Total"]) / thread_scale) * 100.0,
                ]
            )


def make_plots(rows: list[tuple[int, dict[str, float]]], results_dir: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to generate plots") from exc

    results_dir.mkdir(parents=True, exist_ok=True)

    baseline_threads = rows[0][0]
    baseline = rows[0][1]
    thread_list = [t for t, _ in rows]
    labels = ["Stage 1", "Stage 2", "Stage 3", "Total"]
    colors = {
        "Stage 1": "#4C72B0",
        "Stage 2": "#DD8452",
        "Stage 3": "#55A868",
        "Total": "#C44E52",
    }

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    for label in labels:
        values = [data[label] for _, data in rows]
        ax.plot(thread_list, values, marker="o", linewidth=2, markersize=7, label=label, color=colors[label])
    ax.set_xlabel("OMP threads")
    ax.set_ylabel("Time (s)")
    ax.set_title("OMP scaling: stage times")
    ax.set_xticks(thread_list)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "omp_stage_times.png", dpi=180)
    plt.close(fig)

    x = list(range(len(thread_list)))
    width = 0.2
    offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]

    fig, ax = plt.subplots(figsize=(9, 5))
    for offset, label in zip(offsets, labels):
        values = [data[label] for _, data in rows]
        xpos = [xi + offset for xi in x]
        ax.bar(xpos, values, width=width, label=label, color=colors[label])

    ax.set_xlabel("OMP threads")
    ax.set_ylabel("Time (s)")
    ax.set_title("OMP scaling: stage-time histogram")
    ax.set_xticks(x)
    ax.set_xticklabels(thread_list)
    ax.grid(axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "omp_stage_histogram.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for label in labels:
        speedup = [baseline[label] / data[label] for _, data in rows]
        ax.plot(thread_list, speedup, marker="o", linewidth=2, markersize=7, label=label, color=colors[label])
    ideal = [tc / baseline_threads for tc in thread_list]
    ax.plot(thread_list, ideal, "--", color="gray", linewidth=1.5, label="Ideal")
    ax.set_xlabel("OMP threads")
    ax.set_ylabel(f"Speedup vs OMP={baseline_threads}")
    ax.set_title("OMP scaling: stage speedup")
    ax.set_xticks(thread_list)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "omp_stage_speedup.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for label in labels:
        efficiency = [
            ((baseline[label] / data[label]) / (thread_count / baseline_threads)) * 100.0
            for thread_count, data in rows
        ]
        ax.plot(thread_list, efficiency, marker="o", linewidth=2, markersize=7, label=label, color=colors[label])
    ax.axhline(100.0, linestyle="--", color="gray", linewidth=1.5, label="Ideal")
    ax.set_xlabel("OMP threads")
    ax.set_ylabel(f"Efficiency vs OMP={baseline_threads} (%)")
    ax.set_title("OMP scaling: stage efficiency")
    ax.set_xticks(thread_list)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "omp_stage_efficiency.png", dpi=180)
    plt.close(fig)


def print_summary(rows: list[tuple[int, dict[str, float]]]) -> None:
    baseline_threads = rows[0][0]
    baseline = rows[0][1]
    print(f"Baseline OMP threads: {baseline_threads}")
    print("OMP  Stage 1    Stage 2    Stage 3      Total   S1 spd   S2 spd   S3 spd  Tot spd")
    print("-----------------------------------------------------------------------------------")
    for threads, data in rows:
        s1 = baseline["Stage 1"] / data["Stage 1"]
        s2 = baseline["Stage 2"] / data["Stage 2"]
        s3 = baseline["Stage 3"] / data["Stage 3"]
        st = baseline["Total"] / data["Total"]
        print(
            f"{threads:>3}  "
            f"{data['Stage 1']:>7.3f}  "
            f"{data['Stage 2']:>9.3f}  "
            f"{data['Stage 3']:>9.3f}  "
            f"{data['Total']:>9.3f}  "
            f"{s1:>7.3f}  "
            f"{s2:>7.3f}  "
            f"{s3:>7.3f}  "
            f"{st:>7.3f}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=None,
        help="Experiment root containing runs/, results/, metadata/, and scripts/",
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=None,
        help="Directory containing omp-*/ subdirectories with bin-generation logs",
    )
    parser.add_argument("--results", type=Path, default=None, help="Output directory for plots and CSV")
    args = parser.parse_args()

    base_dir = args.base_dir or resolve_default_base_dir()
    runs_dir = args.runs_dir or ((base_dir / "runs") if (base_dir / "runs").is_dir() else base_dir)
    results_dir = args.results or (base_dir / "results")

    logs = discover_logs(runs_dir)
    if not logs:
        print(f"No logs found under {runs_dir}", file=sys.stderr)
        return 1

    rows = collect_rows(logs)
    if not rows:
        print(f"No complete timing summaries found under {runs_dir}", file=sys.stderr)
        return 1

    print_summary(rows)
    results_dir.mkdir(parents=True, exist_ok=True)
    write_csv(rows, results_dir / "omp_stage_summary.csv")
    make_plots(rows, results_dir)

    print(f"\nSaved plots and CSV to {results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
