#!/bin/bash
# Login-node half of the custom multistep build: a hashrun.sh-style checkout of the production
# abacus commit with subprojects/zeldovich_mpi replaced by this repo's working tree (branch
# neutrinos, including uncommitted changes, recorded in provenance/). Needs GitHub access for
# the submodules, so run it on a login node; then qsub pbs_build_multistep_nu.sh (compute node).
#
# Usage: bash prepare_checkout.sh > prepare_checkout.log 2>&1
set -euo pipefail

CODE_HASH=9575f04bdaf2f7c88c2b160ac0d60b428a80b4d3        # production build of the AbacusAurora runs
CODE_REPO=/home/helenshao/InitialConditions/abacus
ZMPI=$CODE_REPO/subprojects/zeldovich_mpi
DEST=${DEST:-/home/helenshao/abacus-store/9575f04_zmpi-neutrinos}

[[ -e $DEST ]] && { echo "$DEST exists; move it aside first" >&2; exit 2; }
mkdir -p "$(dirname "$DEST")"

git init -q "$DEST"
git -C "$DEST" fetch --quiet --depth 1 "$CODE_REPO" "$CODE_HASH"
git -C "$DEST" checkout --quiet --detach FETCH_HEAD
git -C "$DEST" submodule update --quiet --init --recursive --depth 1
echo "abacus: $(git -C "$DEST" rev-parse HEAD)"
echo "production zeldovich_mpi submodule: $(git -C "$DEST/subprojects/zeldovich_mpi" rev-parse HEAD)"

# swap in the neutrinos working tree (tracked files only; the pinned submodule is kept aside)
mv "$DEST/subprojects/zeldovich_mpi" "$DEST/subprojects/zeldovich_mpi.production"
mkdir "$DEST/subprojects/zeldovich_mpi"
(cd "$ZMPI" && git ls-files -z | rsync -a --from0 --files-from=- ./ "$DEST/subprojects/zeldovich_mpi/")

mkdir -p "$DEST/provenance"
{
    echo "date: $(date)"
    echo "abacus: $CODE_HASH"
    echo "zeldovich_mpi: $(git -C "$ZMPI" rev-parse HEAD) ($(git -C "$ZMPI" branch --show-current)) + uncommitted diff below"
    git -C "$ZMPI" status --short -- src deps meson.build meson.options
} > "$DEST/provenance/zeldovich_mpi.txt"
git -C "$ZMPI" diff HEAD -- src deps meson.build meson.options > "$DEST/provenance/zeldovich_mpi_uncommitted.diff"

# everything that needs the network, so the compute-node build is offline: python deps (incl.
# meson, abacusutils from git) without building abacus itself, and the meson wrap tarballs
# (the production store's packagecache is not readable)
(cd "$DEST" && . ./env/aurora.sh && uv sync --no-install-project && meson subprojects download)
ls "$DEST/subprojects" "$DEST/subprojects/packagecache"
echo "prepared $DEST"
