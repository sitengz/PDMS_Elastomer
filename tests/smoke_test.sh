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

run_case base --strand-length 16 --strand-count 12 \
    --crosslinker-length 16
run_case full --filler-wt 5 --filler-length 8 \
    --strand-length 24 --strand-count 12 --crosslinker-length 16 \
    --moderator-count 2 \
    --thickness 40
run_case ratio --strand-length 16 --strand-count 12 \
    --crosslinker-length 16 --stoichiometry 2:1
run_case config --config "$repo_root/tests/basic_model.conf" --strand-count 16
run_case ring_bifunctional --strand-topology ring \
    --strand-length 16 --strand-count 12 --strand-functionality 2 \
    --strand-reactive-distribution regular --crosslinker-length 16 \
    --output data.ring_bifunctional
run_case ring_tetrafunctional --strand-topology ring \
    --strand-length 16 --strand-count 12 --strand-functionality 4 \
    --strand-reactive-distribution random --strand-reactive-seed 12345 \
    --crosslinker-length 16 --output data.ring_tetrafunctional
run_case star_3arm --strand-topology star --strand-length 4 \
    --strand-count 4 --strand-arm-count 3 --crosslinker-length 8 \
    --density 0.02 --output data.star_3arm
run_case star_4arm --strand-topology star --strand-length 4 \
    --strand-count 4 --strand-arm-count 4 --crosslinker-length 8 \
    --density 0.02 --output data.star_4arm
run_case star_6arm --strand-topology star --strand-length 4 \
    --strand-count 2 --strand-arm-count 6 --crosslinker-length 8 \
    --density 0.02 --output data.star_6arm
run_case star_8arm --strand-topology star --strand-length 4 \
    --strand-count 2 --strand-arm-count 8 --crosslinker-length 8 \
    --density 0.02 --output data.star_8arm

base_dir="$test_root/base/PDMS_elastomer"
full_dir="$test_root/full/PDMS_elastomer_filler_N8_5wt_film_Lz40"
ratio_dir="$test_root/ratio/PDMS_elastomer"
config_dir="$test_root/config/from_config"
ring_bifunctional_dir="$test_root/ring_bifunctional/ring_bifunctional"
ring_tetrafunctional_dir="$test_root/ring_tetrafunctional/ring_tetrafunctional"
star_3arm_dir="$test_root/star_3arm/star_3arm"
star_4arm_dir="$test_root/star_4arm/star_4arm"
star_6arm_dir="$test_root/star_6arm/star_6arm"
star_8arm_dir="$test_root/star_8arm/star_8arm"

test -f "$base_dir/data.PDMS_elastomer"
test -f "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_Lz40"
test -f "$ratio_dir/data.PDMS_elastomer"
test -f "$config_dir/data.from_config"
test -f "$ring_bifunctional_dir/data.ring_bifunctional"
test -f "$ring_tetrafunctional_dir/data.ring_tetrafunctional"
test -f "$star_3arm_dir/data.star_3arm"
test -f "$star_4arm_dir/data.star_4arm"
test -f "$star_6arm_dir/data.star_6arm"
test -f "$star_8arm_dir/data.star_8arm"

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

check_data "$base_dir/data.PDMS_elastomer"
check_data "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_Lz40"
check_data "$ratio_dir/data.PDMS_elastomer"
check_data "$config_dir/data.from_config"
check_data "$ring_bifunctional_dir/data.ring_bifunctional"
check_data "$ring_tetrafunctional_dir/data.ring_tetrafunctional"
check_data "$star_3arm_dir/data.star_3arm"
check_data "$star_4arm_dir/data.star_4arm"
check_data "$star_6arm_dir/data.star_6arm"
check_data "$star_8arm_dir/data.star_8arm"

grep -q 'prob 0.10000000' "$base_dir/in.PDMS_elastomer"
test "$(grep -c 'wall/lj126' "$full_dir/in.PDMS_elastomer_filler_N8_5wt_film_Lz40")" -eq 4
grep -q '"format": "pdms-elastomer-model-info"' "$base_dir/PDMS_elastomer.info"
grep -q '"format_version": 3' "$base_dir/PDMS_elastomer.info"
grep -q '"strands": {"component": 1, "N": 16, "M": 12' "$base_dir/PDMS_elastomer.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 6' "$base_dir/PDMS_elastomer.info"
grep -q '"moderators": {"component": 3, "N": 5, "M": 2' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_Lz40.info"
grep -q '"filler": {"component": 4, "N": 8, "M": 3' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_Lz40.info"
grep -q '"molecule_id_start": 19, "molecule_id_end": 20' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_Lz40.info"
grep -q '"molecule_id_start": 21, "molecule_id_end": 23' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_Lz40.info"
grep -q '"moderators_included": false' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_Lz40.info"
grep -q '"extra_moderator_functional_groups": 8' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_Lz40.info"
grep -q '"requested": "2:1"' "$ratio_dir/PDMS_elastomer.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 3' "$ratio_dir/PDMS_elastomer.info"
grep -q '"expected_frames": 1001' "$base_dir/PDMS_elastomer.info"
grep -q '"maximum_allowed_beads": 150000' "$base_dir/PDMS_elastomer.info"
grep -q '"config_file": ' "$config_dir/from_config.info"
grep -q '"strands": {"component": 1, "N": 16, "M": 16' "$config_dir/from_config.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 8' "$config_dir/from_config.info"
grep -q '"topology": "ring"' "$ring_bifunctional_dir/ring_bifunctional.info"
grep -q '"functionality": 2' "$ring_bifunctional_dir/ring_bifunctional.info"
grep -q '"reactive_distribution": "regular"' "$ring_bifunctional_dir/ring_bifunctional.info"
grep -q '"reactive_bead_sites": \[1, 9\]' "$ring_bifunctional_dir/ring_bifunctional.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 6' "$ring_bifunctional_dir/ring_bifunctional.info"
grep -q '^282 bonds$' "$ring_bifunctional_dir/data.ring_bifunctional"
grep -q '"topology": "ring"' "$ring_tetrafunctional_dir/ring_tetrafunctional.info"
grep -q '"functionality": 4' "$ring_tetrafunctional_dir/ring_tetrafunctional.info"
grep -q '"reactive_distribution": "random"' "$ring_tetrafunctional_dir/ring_tetrafunctional.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 12' "$ring_tetrafunctional_dir/ring_tetrafunctional.info"
grep -q '^372 bonds$' "$ring_tetrafunctional_dir/data.ring_tetrafunctional"

awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { in_atoms = 0 }
    in_atoms && NF == 7 && $2 <= 12 && $3 == 2 { ++sites[$2] }
    END {
        for (molecule = 1; molecule <= 12; ++molecule)
            if (sites[molecule] != 4) exit 1
    }
' "$ring_tetrafunctional_dir/data.ring_tetrafunctional"

awk '
    /^Bonds$/ { in_bonds = 1; next }
    in_bonds && NF == 4 && (($3 == 16 && $4 == 1) || ($3 == 1 && $4 == 16)) {
        closure = 1
    }
    END { exit !closure }
' "$ring_bifunctional_dir/data.ring_bifunctional"

check_star() {
    local directory=$1
    local case_name=$2
    local arms=$3
    local centers=$4
    local strand_beads=$5
    local strand_count=$6
    local crosslinker_count=$7
    local center_degree=$8
    local data="$directory/data.$case_name"
    local info="$directory/$case_name.info"

    grep -q '"topology": "star"' "$info"
    grep -q "\"functionality\": $arms" "$info"
    grep -q "\"beads_per_molecule\": $strand_beads" "$info"
    grep -q "\"star_arm_count\": $arms" "$info"
    grep -q "\"star_center_count\": $centers" "$info"
    grep -q '"star_arm_length": 4' "$info"
    grep -q '"reactive_distribution": "arm_ends"' "$info"
    grep -q "\"strands\": {\"component\": 1, \"N\": $strand_beads, \"M\": $strand_count" "$info"
    grep -q "\"crosslinkers\": {\"component\": 2, \"N\": 8, \"M\": $crosslinker_count" "$info"

    awk -v molecules="$strand_count" -v arms="$arms" '
        /^Atoms # full$/ { in_atoms = 1; next }
        in_atoms && /^Bonds$/ { in_atoms = 0 }
        in_atoms && NF == 7 && $2 <= molecules && $3 == 2 {
            ++reactive[$2]
        }
        END {
            for (molecule = 1; molecule <= molecules; ++molecule)
                if (reactive[molecule] != arms) exit 1
        }
    ' "$data"

    awk -v centers="$centers" -v expected="$center_degree" '
        /^Bonds$/ { in_bonds = 1; next }
        in_bonds && /^[A-Za-z]/ { in_bonds = 0 }
        in_bonds && NF == 4 {
            if ($3 <= centers) ++degree[$3]
            if ($4 <= centers) ++degree[$4]
        }
        END {
            for (atom = 1; atom <= centers; ++atom)
                if (degree[atom] != expected) exit 1
        }
    ' "$data"
}

check_star "$star_3arm_dir" star_3arm 3 1 13 4 3 3
check_star "$star_4arm_dir" star_4arm 4 1 17 4 4 4
check_star "$star_6arm_dir" star_6arm 6 2 26 2 3 4
check_star "$star_8arm_dir" star_8arm 8 3 35 2 4 4

limit_dir="$test_root/limit"
mkdir -p "$limit_dir"
if (
    cd "$limit_dir"
    "$generator" --strand-count 1100 >/dev/null 2>limit.err
); then
    echo "oversized model unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'maximum allowed is 150000' "$limit_dir/limit.err"

invalid_ring_dir="$test_root/invalid_ring"
mkdir -p "$invalid_ring_dir"
if (
    cd "$invalid_ring_dir"
    "$generator" --strand-topology ring --strand-length 15 \
        --strand-functionality 2 --strand-reactive-distribution regular \
        >/dev/null 2>ring.err
); then
    echo "invalid regular ring unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'divisible by functionality' "$invalid_ring_dir/ring.err"

invalid_star_dir="$test_root/invalid_star"
mkdir -p "$invalid_star_dir"
if (
    cd "$invalid_star_dir"
    "$generator" --strand-topology star --strand-arm-count 5 \
        >/dev/null 2>star.err
); then
    echo "unsupported star arm count unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'Star arm count must be 3, 4, 6, or 8' "$invalid_star_dir/star.err"

fixed_geometry_dir="$test_root/fixed_geometry"
mkdir -p "$fixed_geometry_dir"
for removed_option in --bond-length --spacing; do
    if (
        cd "$fixed_geometry_dir"
        "$generator" "$removed_option" 3 >/dev/null 2>option.err
    ); then
        echo "$removed_option unexpectedly remained configurable" >&2
        exit 1
    fi
    grep -q "Unknown option: $removed_option" "$fixed_geometry_dir/option.err"
done

awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { in_atoms = 0 }
    in_atoms && NF == 7 {
        molecule = $2
        type = $3
        if (molecule <= 12 && type == 2) ++strand_sites[molecule]
        else if (molecule >= 13 && molecule <= 18 && type == 3) ++crosslinker_sites[molecule]
        else if (molecule >= 19 && molecule <= 20 && type == 2) ++moderator_sites[molecule]
        else if (molecule >= 21 && type != 1) bad_filler = 1
    }
    END {
        for (molecule = 1; molecule <= 12; ++molecule)
            if (strand_sites[molecule] != 2) exit 1
        for (molecule = 13; molecule <= 18; ++molecule)
            if (crosslinker_sites[molecule] != 4) exit 1
        for (molecule = 19; molecule <= 20; ++molecule)
            if (moderator_sites[molecule] != 4) exit 1
        if (bad_filler) exit 1
    }
' "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_Lz40"

echo "PDMS elastomer smoke tests passed"
