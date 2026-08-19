#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
model_generator="$repo_root/bin/pdms_elastomer_generator"
tensile_generator="$repo_root/bin/tensile_test_generator"
test_root=$(mktemp -d)
trap 'rm -rf -- "$test_root"' EXIT

generate_source() {
    local name=$1
    shift
    (
        cd "$test_root"
        "$model_generator" --strand-length 12 --strand-count 8 \
            --crosslinker-length 8 --output "data.$name" "$@" >/dev/null
    )
    cp "$test_root/$name/data.$name" "$test_root/$name/data.$name.npt_eq"
}

generate_source tensile_bulk
generate_source tensile_film --thickness 40

"$tensile_generator" \
    "$test_root/tensile_bulk/data.tensile_bulk.npt_eq" \
    "$test_root/tensile_bulk/tensile_bulk.info" >/dev/null
"$tensile_generator" \
    "$test_root/tensile_film/data.tensile_film.npt_eq" \
    "$test_root/tensile_film/tensile_film.info" >/dev/null

for direction in x y z; do
    test -f "$test_root/tensile_bulk/tensile/$direction/in.tensile.tensile_bulk.$direction"
    test -f "$test_root/tensile_bulk/tensile/$direction/tensile.tensile_bulk.$direction.info"
    test -x "$test_root/tensile_bulk/tensile/$direction/submit.tensile.tensile_bulk.$direction.sh"
done

for direction in x y; do
    input="$test_root/tensile_film/tensile/$direction/in.tensile.tensile_film.$direction"
    test -f "$input"
    grep -q '^boundary        p p f$' "$input"
    test "$(grep -c 'wall/lj126' "$input")" -eq 2
    grep -q '^comm_modify     cutoff 25$' "$input"
    grep -q '^pair_style      lj/gromacs 12 15$' "$input"
    grep -q '^run             5000000$' "$input"
    grep -q '^run             20000000$' "$input"
    grep -q "^fix             deform_box all deform 1 $direction erate 2.000000000000e-09" "$input"
    grep -q '^fix             samples all ave/time 100 10 1000' "$input"
    if grep -Eq '^(minimize|dump|fix +xlink)' "$input"; then
        echo "Tensile input unexpectedly minimizes, dumps, or crosslinks" >&2
        exit 1
    fi
done

test ! -e "$test_root/tensile_film/tensile/z"
grep -q 'variable        reduced_strain equal v_lambda_x\^2-v_lambda_y\^2' \
    "$test_root/tensile_film/tensile/x/in.tensile.tensile_film.x"
grep -q '"stress_volume_correction": "box volume / material volume = Lz/H"' \
    "$test_root/tensile_film/tensile/x/tensile.tensile_film.x.info"
grep -q '"default_fit_engineering_strain_range": \[0.0, 0.05\]' \
    "$test_root/tensile_film/tensile/x/tensile.tensile_film.x.info"

override_output="$test_root/wrapper_output"
NPT_EQ_FILE="$test_root/tensile_film/data.tensile_film.npt_eq" \
INFO_FILE="$test_root/tensile_film/tensile_film.info" \
TENSILE_OUTPUT_DIR="$override_output" \
    bash "$repo_root/simulations/01_linear_reference/film_2Ree/run_tensile.sh" \
        --direction x >/dev/null
test -f "$override_output/x/in.tensile.tensile_film.x"

echo "Tensile generator tests passed"
