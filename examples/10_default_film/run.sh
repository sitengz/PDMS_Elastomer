#!/usr/bin/env bash
set -euo pipefail

example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$example_dir/../.." && pwd)
film_thickness=${FILM_THICKNESS_ANGSTROM:-}

if [[ -z "$film_thickness" ]]; then
    echo "Set FILM_THICKNESS_ANGSTROM to the desired nominal wall-free film thickness." >&2
    echo "Example: FILM_THICKNESS_ANGSTROM=100 bash examples/10_default_film/run.sh" >&2
    exit 2
fi

make -C "$repo_root" bin/pdms_elastomer_generator
cd "$example_dir"
"$repo_root/bin/pdms_elastomer_generator" \
    --config model.conf \
    --thickness "$film_thickness" \
    "$@"
