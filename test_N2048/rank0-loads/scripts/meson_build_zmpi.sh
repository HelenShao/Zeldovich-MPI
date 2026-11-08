# Source after: module load … (same as PBS).
# Sets: ZMPI_EXE, WISDOM_EXE (Meson output paths).
# Repo root: directory containing top-level meson.build — either this script’s dir
# (if placed in repo root) or three levels up from test_N2048/rank0-loads/scripts/.
_SCRIPT_PATH="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd "$(dirname "$_SCRIPT_PATH")" && pwd)"

if [[ -n "${REPO_ROOT:-}" ]]; then
  :
elif [[ -f "$_SCRIPT_DIR/meson.build" ]]; then
  REPO_ROOT="$_SCRIPT_DIR"
elif [[ -f "$_SCRIPT_DIR/../../../meson.build" ]]; then
  REPO_ROOT="$(cd "$_SCRIPT_DIR/../../.." && pwd)"
else
  echo "meson_build_zmpi.sh: cannot find repo root (meson.build). Set REPO_ROOT." >&2
  exit 1
fi

MESON_BUILD="${MESON_BUILD:-$REPO_ROOT/build}"

cd "$REPO_ROOT" || { echo "meson_build_zmpi.sh: cd $REPO_ROOT failed"; exit 1; }

MESON_SETUP=(
  -Dproduction_mode=true
  -Dparticle_output_mode=4
  -Dskip_file_write=0
  -Dparallelize_z_loop=1
  -Ddebug_prints=false
  -Ddebug_rng_consistency=false
  -Ddebug_rng_skip=false
  -Duse_fftw_wisdom=true
  -Dextra_cpp_args=-g
)

if [[ -f "$MESON_BUILD/meson-private/coredata.dat" ]]; then
  meson setup "$MESON_BUILD" "${MESON_SETUP[@]}" --reconfigure
else
  meson setup "$MESON_BUILD" "${MESON_SETUP[@]}"
fi
meson compile -C "$MESON_BUILD"

export ZMPI_EXE="$MESON_BUILD/src/Zeldovich_MPI"
export WISDOM_EXE="$MESON_BUILD/src/wisdom_rank0"
if [[ ! -x "$ZMPI_EXE" ]]; then
  echo "Meson build failed: missing or non-executable $ZMPI_EXE"
  exit 1
fi
if [[ ! -x "$WISDOM_EXE" ]]; then
  echo "Meson build failed: missing or non-executable $WISDOM_EXE"
  exit 1
fi
