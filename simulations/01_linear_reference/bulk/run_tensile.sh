#!/usr/bin/env bash
set -euo pipefail

geometry_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "$geometry_dir/../../run_tensile_case.sh" "$geometry_dir" "$@"
