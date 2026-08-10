#!/usr/bin/env bash
set -euo pipefail

example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$example_dir/../.." && pwd)

make -C "$repo_root" bin/pdms_elastomer_generator
cd "$example_dir"
"$repo_root/bin/pdms_elastomer_generator" --config model.conf "$@"
