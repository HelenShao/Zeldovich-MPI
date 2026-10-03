#!/bin/bash
#PBS -N nu_st5_emb
#PBS -A Abacus
#PBS -q debug
#PBS -l select=2
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/stage5/pbs_stage5_embedded.log

# Stage 5, embedded half (2 nodes = 4 ranks: NumZRanks = 2 needs >= 2 X ranks): the custom multistep runs the real FLAMINGO_nu_c202.par2 (only NP, CPD and the
# IC-step controls overridden), generating ICs from the compiled-in FLAMINGO noise directory. Then:
#   - the "[external noise] reading" log line must name the FLAMINGO noise
#   - the ZD_* keys Abacus wrote must equal the hand-written standalone par
#   - ICs vs pbs_stage5_standalone.sh (FFT-roundoff tolerance: the builds plan FFTs differently)

set -o pipefail
HERE=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/stage5
CMP=$HERE/../roundtrip256/compare_rt256.py
PAR2=/home/helenshao/InitialConditions/AbacusAurora/Simulations/FLAMINGO_nu_c202.par2
CODE=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos
BASE=/flare/Abacus/helenshao/stage5_c202_N300
SA=$BASE/standalone
EMB=$BASE/embedded
FFT_TOL=${FFT_TOL:-1e-5}

if [[ ! -e /opt/aurora/26.181.0/oneapi/compiler/latest/lib/libsycl.so.9 ]]; then
    echo "wrong system image on $(hostname): /opt/aurora/26.181.0 missing (needs the 26.181.0 compute image)" >&2
    exit 3
fi
[[ -x $CODE/build/src/Multistep/multistep ]] || { echo "custom multistep missing (custom_build job)" >&2; exit 1; }
[[ -d $SA/ic ]] || { echo "no standalone ICs in $SA (pbs_stage5_standalone.sh)" >&2; exit 1; }
[[ -e $EMB ]] && { echo "$EMB exists; move it aside first" >&2; exit 2; }
mkdir -p "$EMB/ic"

(
    cd "$CODE" && . ./env/aurora.sh
    export ABACUS_WORKING_ROOT=$EMB ABACUS_OUTPUT_ROOT=$EMB ABACUS_CHECKPOINT_ROOT=$EMB
    unset ZD_EXTERNAL_NOISE_DIR
    {
        date; hostname
        echo "code: $CODE"; ls "$CODE/provenance" 2>/dev/null
        ls -l "$ABACUS_BUILD/src/Multistep/multistep"
        echo "par2: $PAR2 (AbacusAurora $(git -C /home/helenshao/InitialConditions/AbacusAurora rev-parse --short HEAD))"
    } > "$EMB/provenance.txt"
    cd "$EMB"
    t0=$SECONDS
    python -u -m abacus.run "$PAR2" -n 1 -P "NP=300**3" -P CPD=125 \
        -P ICCleanupMode=Never -P InitialConditionsDirectory="$EMB/ic" -P UseSCR=0 > run.log 2>&1
    echo "abacus.run rc=$? after $((SECONDS - t0)) s"
)

{
    echo "Stage 5: embedded (custom multistep, FLAMINGO_nu_c202.par2) vs standalone, 300^3 ($(date))"
    echo
    n=$(grep -rh "\[external noise\] reading" "$EMB" --include='*.log' 2>/dev/null | head -1)
    echo "[log] embedded: ${n:-<none>}"
    [[ $n == *flamingo_panphasia_N2048/noise* ]] && echo "[log] -> PASS" || echo "[log] -> FAIL"
    echo
    echo "[par] ZD_* / IC keys, hand-written standalone par vs Abacus's abacus.par:"
    python3 - "$SA/param_stage5_c202_N300.par" "$(find "$EMB" -name abacus.par | head -1)" <<'EOF'
import re, sys
KEYS = ("BoxSize", "NP", "CPD", "NumZRanks", "InitialRedshift", "Omega_M", "ICFormat")
def load(fn):
    d = {}
    for line in open(fn):
        m = re.match(r"^\s*(\w+)\s*=\s*(.*?)\s*$", line.split("#")[0])
        if m and (m[1].startswith("ZD_") or m[1] in KEYS):
            d[m[1]] = m[2]
    return d
def same(x, y):
    if x == y:
        return True
    try:
        return abs(float(eval(x)) - float(eval(y))) <= 1e-12 * max(1.0, abs(float(eval(x))))
    except Exception:
        return False
a, b = load(sys.argv[1]), load(sys.argv[2])
bad = []
for k in sorted(set(a) | set(b)):
    if k == "ZD_PLT_filename":
        continue
    if k not in a or k not in b:
        print(f"  only in {'embedded' if k in b else 'standalone'}: {k} = {b.get(k, a.get(k))}")
        continue
    if not same(a[k], b[k]):
        bad.append(k)
        print(f"  differs: {k}: standalone={a[k]} embedded={b[k]}")
import filecmp
pa, pb = a["ZD_PLT_filename"].strip('"'), b.get("ZD_PLT_filename", '""').strip('"')
plt_same = filecmp.cmp(pa, pb, shallow=False) if pa != pb else True
print(f"  ZD_PLT_filename: {pb} ({'same content as standalone' if plt_same else 'DIFFERENT eigmodes file'})")
print(f"  {len(set(a) & set(b))} common keys, {len(bad)} differ -> {'PASS' if not bad and plt_same else 'FAIL'}")
EOF
    echo
} > "$HERE/results.txt" 2>&1

python3 "$CMP" --N 300 --noise-dir "" --ic-plain "" --d-ref "" --d-test "" \
    --ic-ref "$SA/ic" --ic-test "$EMB/ic" --fft-tol "$FFT_TOL" \
    --out "$BASE/compare_results.txt" > "$BASE/compare.log" 2>&1
rc=$?
cat "$BASE/compare_results.txt" >> "$HERE/results.txt"
echo "compare rc=$rc" >> "$HERE/results.txt"
cat "$HERE/results.txt"
exit $rc
