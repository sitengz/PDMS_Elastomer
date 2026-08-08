#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
generator="$repo_root/bin/pdms_elastomer_generator"
test_root=$(mktemp -d)
trap 'rm -rf -- "$test_root"' EXIT

run_case() {
    local name=$1
    shift
    mkdir -p "$test_root/$name"
    (
        cd "$test_root/$name"
        "$generator" "$@" >/dev/null 2>&1
    )
}

run_case v22 --formulation V22 --n1 16 --m1 12 --n2 16 \
    --functionality 8 --spacing 3
run_case v35 --filler-wt 5 --formulation V35 --filler-length 8 \
    --n1 24 --m1 12 --n2 16 --functionality 4 --spacing 3 \
    --thickness 40

v22_dir="$test_root/v22/V22_no_filler"
v35_dir="$test_root/v35/V35_PDMS_N8_5wt_film_Lz40"

test -f "$v22_dir/data.V22_no_filler"
test -f "$v35_dir/data.V35_PDMS_N8_5wt_film_Lz40"

check_data() {
    local data=$1
    awk '
        /^3 atom types$/ { atoms = 1 }
        /^2 bond types$/ { bonds = 1 }
        /^Bonds$/ { in_bonds = 1; next }
        in_bonds && /^[A-Za-z]/ { in_bonds = 0 }
        in_bonds && NF == 4 && $2 != 1 { bad_bond = 1 }
        END { exit !(atoms && bonds && !bad_bond) }
    ' "$data"
}

check_data "$v22_dir/data.V22_no_filler"
check_data "$v35_dir/data.V35_PDMS_N8_5wt_film_Lz40"

grep -q 'prob 0.10000000' "$v22_dir/in.V22_no_filler"
grep -q 'prob 0.50000000' "$v35_dir/in.V35_PDMS_N8_5wt_film_Lz40"
test "$(grep -c 'wall/lj126' "$v35_dir/in.V35_PDMS_N8_5wt_film_Lz40")" -eq 4
grep -q '"expected_frames": 1001' "$v22_dir/V22_no_filler.info"

"$repo_root/bin/binary_morphology_analyzer" \
    "$v35_dir/data.V35_PDMS_N8_5wt_film_Lz40" \
    "$v35_dir/V35_PDMS_N8_5wt_film_Lz40.info" \
    --grid-spacing 8 --contact-cutoff 8 \
    --output "$test_root/morphology_metrics.dat" \
    --structure-output "$test_root/morphology_structure.dat" \
    --field-output "$test_root/morphology_field.dat" \
    --report-output "$test_root/morphology_report.txt" >/dev/null
awk 'NF != 22 { exit 1 }' "$test_root/morphology_metrics.dat"
awk 'NF != 10 { exit 1 }' "$test_root/morphology_field.dat"
test -s "$test_root/morphology_structure.dat"
grep -q 'BINARY ELASTOMER/FILLER MORPHOLOGY ANALYSIS' \
    "$test_root/morphology_report.txt"

echo "PDMS elastomer smoke tests passed"
