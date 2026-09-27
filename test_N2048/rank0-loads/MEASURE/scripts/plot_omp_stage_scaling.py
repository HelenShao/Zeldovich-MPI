#!/usr/bin/env python3
"""Delegate to the shared plotter under rank0-loads/scripts/."""

from pathlib import Path
import runpy

_src = Path(__file__).resolve().parents[2] / "scripts" / "plot_omp_stage_scaling.py"
if not _src.is_file():
    raise FileNotFoundError(f"Expected plotter at {_src}")
runpy.run_path(str(_src), run_name="__main__")
