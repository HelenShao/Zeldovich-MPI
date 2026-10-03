#!/bin/bash
#PBS -N nu256_smoke
#PBS -A Abacus
#PBS -q capacity
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/smoke/pbs_smoke_run.log

# Abacus smoke test (ABACUS_SMOKE_TEST.md stage 3): the production multistep reads the
# external-noise ICs (ExternalICs = 1) and runs 256^3 cosm202 from z = 31 to 0.

set -o pipefail

CODE=/home/eisenste/abacus-store/9575f04bdaf2f7c88c2b160ac0d60b428a80b4d3   # production build (AbacusAurora emulator runs)
PAR2=/home/helenshao/InitialConditions/AbacusAurora/Simulations/smoke_nu256_c202.par2
SMOKE=/home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/smoke
ICD=/flare/Abacus/helenshao/ICs/nu256_c202_panphasia
OUT=/flare/Abacus/helenshao/AbacusAurora_smoke_nu256_c202

cd "$CODE" && . ./env/aurora.sh
cd "$SMOKE"

# the production binary links the 26.181.0 stack (libsycl.so.9); the default image only has .so.8
if [[ ! -e /opt/aurora/26.181.0/oneapi/compiler/latest/lib/libsycl.so.9 ]]; then
    echo "wrong system image on $(hostname): /opt/aurora/26.181.0 missing (needs the 26.181.0 compute image)" >&2
    exit 3
fi

if [[ -e $OUT ]]; then
    echo "$OUT already exists; move it aside before a fresh run" >&2
    exit 2
fi

{
    date
    hostname
    echo "code: $CODE"
    echo "par2: $PAR2 (AbacusAurora $(git -C /home/helenshao/InitialConditions/AbacusAurora rev-parse --short HEAD))"
    echo "zeldovich_mpi (IC generator): $(git -C /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi rev-parse HEAD)"
    module list 2>&1
} > provenance.txt

(cd "$ICD/ic" && sha256sum -c ../ic_sha256_source.txt) > ic_sha256_pre.txt 2>&1 \
    || { echo "IC checksums failed before the run" >&2; exit 1; }

echo "=== abacus.run start: $(date) ==="
t0=$SECONDS
# LPT=1 (qsub -v LPT=1): variant S1b. These ICs have vel = displ exactly (f_cluster = 1), so the 2LPT
# velocity bound lpt_vel_scale is pure double roundoff and the float32 kick overflows its 16-bit packing
# (particlestruct.cpp _pack_float assertion); production ICs (f_cluster < 1) are unaffected.
echo "LagrangianPTOrder override: ${LPT:-none (par2 value)}" >> provenance.txt
python -u -m abacus.run "$PAR2" ${LPT:+-P LagrangianPTOrder=$LPT} > run.log 2>&1
rc=$?
echo "=== abacus.run rc=$rc after $((SECONDS - t0)) s: $(date) ==="

(cd "$ICD/ic" && sha256sum -c ../ic_sha256_source.txt) > ic_sha256_post.txt 2>&1
echo "IC checksums after run: $(grep -c ': OK$' ic_sha256_post.txt) OK of $(wc -l < "$ICD/ic_sha256_source.txt")"

grep -h "Embedded Zeldovich_MPI disabled" "$OUT"/log/*.log 2>/dev/null | head -2
ls -d "$OUT"/slice* 2>/dev/null
exit $rc
