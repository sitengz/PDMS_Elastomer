#!/usr/bin/env bash
set -euo pipefail

geometry_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$geometry_dir/../../.." && pwd)

make -C "$repo_root" bin/pdms_elastomer_generator
cd "$geometry_dir"
"$repo_root/bin/pdms_elastomer_generator" \
    --config model.conf \
    "$@"
