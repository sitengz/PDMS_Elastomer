#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 GEOMETRY_DIRECTORY [layer-dynamics-generator options]" >&2
    exit 2
fi

geometry_dir=$(cd -- "$1" && pwd)
shift
repo_root=$(cd -- "$geometry_dir/../../.." && pwd)
model_config="$geometry_dir/model.conf"

if [[ ! -f "$model_config" ]]; then
    echo "Missing model configuration: $model_config" >&2
    exit 2
fi

data_basename=$(awk -F= '
    /^[[:space:]]*output[[:space:]]*=/ {
        value = $2
        sub(/^[[:space:]]+/, "", value)
        sub(/[[:space:]]+$/, "", value)
        print value
    }
' "$model_config")

if [[ -z "$data_basename" ]]; then
    echo "Could not determine output from $model_config" >&2
    exit 2
fi

case_name=${data_basename#data.}
default_case_dir="$geometry_dir/$case_name"
npt_eq=${NPT_EQ_FILE:-"$default_case_dir/data.$case_name.npt_eq"}
model_info=${INFO_FILE:-"$default_case_dir/$case_name.info"}
output_dir=${LAYER_DYNAMICS_OUTPUT_DIR:-"$(dirname -- "$npt_eq")/layer_dynamics"}

if [[ ! -f "$npt_eq" ]]; then
    echo "Missing final equilibrated network: $npt_eq" >&2
    echo "Run and complete the original simulation before run_layer_dynamics.sh." >&2
    echo "For a nonstandard location, set NPT_EQ_FILE explicitly." >&2
    exit 2
fi

if [[ ! -f "$model_info" ]]; then
    echo "Missing matching model information: $model_info" >&2
    echo "For a nonstandard location, set INFO_FILE explicitly." >&2
    exit 2
fi

make -C "$repo_root" bin/layer_dynamics_generator

"$repo_root/bin/layer_dynamics_generator" \
    "$npt_eq" \
    "$model_info" \
    --output-dir "$output_dir" \
    "$@"

echo
echo "Layer-dynamics files are ready in: $output_dir"
find "$output_dir" -maxdepth 1 -name 'submit.layer_dynamics.*.sh' -print | sort
