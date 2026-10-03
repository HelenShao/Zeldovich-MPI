#!/bin/bash
#PBS -N nu_rv2_emb
#PBS -A Abacus
#PBS -q debug-scaling
#PBS -l select=4
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/roundtrip_rockv2/pbs_embedded.log

# rockV2 round trip, part 2 (background, non-gating). Same node count as the production run (4 nodes,
# 2 ranks per node, NumZRanks = 2). Both runs use the unchanged rockV2 par2 and stop after the IC step:
#   ref:  production multistep (eisenste 9575f04), seed noise, embedded Zeldovich-MPI
#   test: custom multistep (9575f04 + zeldovich_mpi neutrinos, LOAD_EXTERNAL_NOISE) reading the
#         seed noise dumped by pbs_seed_dump.sh through ZD_EXTERNAL_NOISE_DIR
# ICs are kept on Flare (ICCleanupMode = Never). Then compare_rt256.py:
#   ref vs test (embedded load round trip) and standalone seed vs ref (standalone == embedded).
# Different binaries use different FFTW plans, so the IC comparison uses an FFT-roundoff tolerance.

set -o pipefail
HERE=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/roundtrip_rockv2
CMP=$HERE/../roundtrip256/compare_rt256.py
PAR2=/home/helenshao/InitialConditions/AbacusAurora/Simulations/rockstar_v2/AbacusAurora_rockV2_c000_ph390_n800.par2
PROD_PAR=/lus/flare/projects/Abacus/AbacusAurora/AbacusAurora_rockV2_c000_ph390_n800/abacus.par
CODE_REF=/home/eisenste/abacus-store/9575f04bdaf2f7c88c2b160ac0d60b428a80b4d3
CODE_TEST=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos
BASE=/flare/Abacus/helenshao/roundtrip_rockv2
SEED=$BASE/seed_standalone
FFT_TOL=${FFT_TOL:-1e-5}

if [[ ! -e /opt/aurora/26.181.0/oneapi/compiler/latest/lib/libsycl.so.9 ]]; then
    echo "wrong system image on $(hostname): /opt/aurora/26.181.0 missing (needs the 26.181.0 compute image)" >&2
    exit 3
fi
[[ -x $CODE_TEST/build/src/Multistep/multistep ]] || { echo "custom multistep missing (custom_build job)" >&2; exit 1; }
[[ -d $SEED/noise_slabs ]] || { echo "no dumped noise in $SEED (pbs_seed_dump.sh)" >&2; exit 1; }
for v in ref test; do
    [[ -e $BASE/$v ]] && { echo "$BASE/$v exists; move it aside first" >&2; exit 2; }
done

run_ic_step() {   # $1 = ref|test, $2 = code dir
    local v=$1 code=$2
    mkdir -p "$BASE/$v/ic"
    (
        cd "$code" && . ./env/aurora.sh
        export ABACUS_WORKING_ROOT=$BASE/$v ABACUS_OUTPUT_ROOT=$BASE/$v ABACUS_CHECKPOINT_ROOT=$BASE/$v
        [[ $v == test ]] && export ZD_EXTERNAL_NOISE_DIR=$SEED/noise_slabs
        {
            date; hostname
            echo "code: $code"; ls "$code/provenance" 2>/dev/null
            ls -l "$ABACUS_BUILD/src/Multistep/multistep"
            echo "par2: $PAR2 (AbacusAurora $(git -C /home/helenshao/InitialConditions/AbacusAurora rev-parse --short HEAD))"
            echo "ZD_EXTERNAL_NOISE_DIR=${ZD_EXTERNAL_NOISE_DIR:-<unset>}"
        } > "$BASE/$v/provenance.txt"
        cd "$BASE/$v"
        python -u -m abacus.run "$PAR2" -n 1 \
            -P ICCleanupMode=Never -P InitialConditionsDirectory="$BASE/$v/ic" \
            -P UseSCR=0 > run.log 2>&1
    )
}

for v in ref test; do
    code=$CODE_REF; [[ $v == test ]] && code=$CODE_TEST
    echo "=== $v start $(date)"
    t0=$SECONDS
    run_ic_step $v "$code"
    echo "=== $v rc=$? after $((SECONDS - t0)) s"
done

{
    echo "rockV2_n800 embedded round trip ($(date))"
    echo
    echo "[log] external-noise lines (test must read $SEED/noise_slabs; ref must have none):"
    for v in ref test; do
        n=$(grep -rh "\[external noise\] reading" "$BASE/$v" --include='*.log' 2>/dev/null | head -1)
        echo "  $v: ${n:-<none>}"
    done
    echo
    echo "[par] ZD_* / IC keys of the ref run vs the finished production abacus.par:"
    python3 - "$PROD_PAR" "$(find "$BASE/ref" -name abacus.par | head -1)" <<'EOF'
import re, sys
def load(fn):
    d = {}
    for line in open(fn):
        m = re.match(r"^\s*(\w+)\s*=\s*(.*?)\s*$", line.split("#")[0])
        if m and (m[1].startswith("ZD_") or m[1] in ("BoxSize", "NP", "CPD", "NumZRanks", "InitialRedshift", "Omega_M", "Omega_Smooth", "ICFormat", "LagrangianPTOrder")):
            d[m[1]] = m[2]
    return d
a, b = load(sys.argv[1]), load(sys.argv[2])
diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
for k in diff:
    print(f"  differs: {k}: prod={a.get(k)} ref={b.get(k)}")
print(f"  {len(a)} keys compared, {len(diff)} differ (paths may differ by store only)")
EOF
} > "$HERE/results.txt" 2>&1

module load frameworks 2>/dev/null
python3 "$CMP" --N 800 --noise-dir "$SEED/noise_slabs" \
    --ic-plain "" --d-ref "" --d-test "" \
    --ic-ref "$BASE/ref/ic" --ic-test "$BASE/test/ic" --fft-tol "$FFT_TOL" \
    --ic-extra "standalone seed vs embedded production" "$SEED/ic" "$BASE/ref/ic" \
    --out "$BASE/compare_results.txt" > "$BASE/compare.log" 2>&1
rc=$?
cat "$BASE/compare_results.txt" >> "$HERE/results.txt"
echo "compare rc=$rc" >> "$HERE/results.txt"
cat "$HERE/results.txt"
exit $rc
