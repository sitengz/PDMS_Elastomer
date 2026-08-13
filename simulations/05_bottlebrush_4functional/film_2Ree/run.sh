#!/usr/bin/env bash
set -euo pipefail

geometry_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$geometry_dir/../../.." && pwd)
bulk_ree=${BULK_REE_ANGSTROM:-}
factor=2

if [[ -z "$bulk_ree" ]]; then
    echo "Set BULK_REE_ANGSTROM to the mean bulk effective-strand Ree in angstrom." >&2
    echo "Example: BULK_REE_ANGSTROM=42.5 bash simulations/05_bottlebrush_4functional/film_2Ree/run.sh" >&2
    exit 2
fi

if ! awk -v value="$bulk_ree" 'BEGIN { exit !(value ~ /^[0-9]+([.][0-9]+)?$/ && value > 0) }'; then
    echo "BULK_REE_ANGSTROM must be a positive number; received: $bulk_ree" >&2
    exit 2
fi

film_thickness=$(awk -v value="$bulk_ree" -v scale="$factor" 'BEGIN { printf "%.8g", value * scale }')

make -C "$repo_root" bin/pdms_elastomer_generator
cd "$geometry_dir"
"$repo_root/bin/pdms_elastomer_generator" \
    --config model.conf \
    --thickness "$film_thickness" \
    "$@"
