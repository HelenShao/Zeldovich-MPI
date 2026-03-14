#!/bin/bash
# Run compare_particle_ics.py for ic_0 through ic_224, saving each to separate PNGs:
#   visualizations/displacement_comparison_ic_0.png, velocity_comparison_ic_0.png, ...
# Run from test_N1024/CPD-225:
#   ./run_compare_particle_ics_all.sh
# Or: bash run_compare_particle_ics_all.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

HERM_BASE="particle_ics"
ZELD_BASE="/home/helenshao/InitialConditions/zeldovich-PLT/output_N1024_CPD_225"
N=1024
OUT_DIR="visualizations"

# Start with a fresh status file so this run overwrites (each ic appends)
rm -f "$OUT_DIR/comparison_status.txt"

for i in $(seq 0 224); do
  echo "===== Comparing ic_$i ====="
  python3 compare_particle_ics.py \
    "${HERM_BASE}/ic_${i}" \
    "${ZELD_BASE}/ic_${i}" \
    "$N" \
    "$OUT_DIR" \
    all \
    "ic_${i}"
done

echo "Done. Check $OUT_DIR for displacement_comparison_ic_*.png, velocity_comparison_ic_*.png, and comparison_status.txt"
