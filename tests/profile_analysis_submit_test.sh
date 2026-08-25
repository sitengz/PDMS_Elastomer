#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
launcher="$repo_root/simulations/submit_profile_analysis.sh"
test_root=$(mktemp -d)
trap 'rm -rf -- "$test_root"' EXIT

bash -n "$launcher"

case_dir="$test_root/demo_case"
layer_dir="$case_dir/layer_dynamics"
mkdir -p "$layer_dir"
touch "$case_dir/data.demo.npt_eq"
touch "$case_dir/demo.info"
touch "$layer_dir/dump.debye_waller.demo.lammpstrj"
touch "$layer_dir/dump.layer_dynamics.demo.lammpstrj"

output=$(bash "$launcher" --dry-run \
    --time-origin-count 20 --time-origin-stride 50 \
    --frame-stride 2 --bin-width 7.5 "$test_root")

grep -q 'DRY RUN: sbatch' <<<"$output"
grep -q -- '--job-name=demo_tamsd' <<<"$output"
grep -q ' 20 50 2 7[.]5 demo$' <<<"$output"
grep -q 'Validated 1 analysis job(s); nothing was submitted.' <<<"$output"
test -d "$case_dir/analysis_demo"

echo "Profile-analysis Slurm launcher tests passed"
