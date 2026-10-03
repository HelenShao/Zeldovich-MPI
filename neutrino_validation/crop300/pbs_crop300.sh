#!/bin/bash
#PBS -N nu_crop300
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/crop300/pbs_crop300.log
#
# Crop test (a): monofonIC GridRes = 300 (Panphasia N0 = 512) and Zeldovich-MPI N = 300, so the
# loader crops 512 -> 300, the production geometry (2048 -> 1800). Then tests 1, 1b, 2, 3, 4.
# Needs build_rt_load (LOAD_EXTERNAL_NOISE + DUMP_D_SLABS, from roundtrip256/build_rt256.sh).
# Refuses to overwrite: move run/ aside for a rerun. Submit: qsub crop300/pbs_crop300.sh
set -eo pipefail

VAL=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation
HERE=$VAL/crop300
ZMPI_EXE=$VAL/../build_rt_load/src/Zeldovich_MPI
MONO_EXE=$HOME/monofonic/build/monofonIC
N=300
N0=512

module load frameworks fftw/3.3.10 hdf5/1.14.6
export NUMEXPR_MAX_THREADS=256
export OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close

[[ -x $ZMPI_EXE ]] || { echo "missing $ZMPI_EXE (run roundtrip256 build first)" >&2; exit 1; }
[[ -e $HERE/run ]] && { echo "$HERE/run exists; move it aside first" >&2; exit 2; }
mkdir -p "$HERE/run/mono" "$HERE/run/external_noise"
cd "$HERE/run"
echo "host $(hostname)  job ${PBS_JOBID:-none}  $(date)"

echo "[1] monofonIC GridRes=$N"
cp "$HERE/monofonic_N300_crop.conf" mono/
(cd mono && mpiexec -n 4 --ppn 4 --depth=8 --cpu-bind depth "$MONO_EXE" monofonic_N300_crop.conf > monofonic.log 2>&1)
mv mono/output_k_space_slab_y.* external_noise/
echo "    slabs: $(ls external_noise/output_k_space_slab_y.* | wc -l) (expect $((N0 / 2 + 1))), size $(stat -c %s external_noise/output_k_space_slab_y.0) B (expect $((N0 * N0 * 16)))"

echo "[2] Zeldovich-MPI N=$N"
cp "$HERE/param_N300_crop.par" .
mpiexec -n 4 --ppn 4 --depth=8 --cpu-bind depth "$ZMPI_EXE" param_N300_crop.par > zmpi_run.log 2>&1
echo "    D_slabs: $(ls D_slabs | wc -l) (expect $((N / 2 + 1)))"

echo "[3] validate_external_noise.py"
python3 "$VAL/validate_external_noise.py" \
    --slab-dir external_noise --N0 $N0 \
    --ours-slabs D_slabs --N $N \
    --box 673.2 --pk "$VAL/pk/cosm202_z31.z1_pk_cb.dat" \
    --mono-hdf5 mono/mono_N300_delta.hdf5 \
    --ic-dir zmpi_out \
    --out validation_results.txt > validate.log 2>&1
cat validation_results.txt
echo "done $(date)"
