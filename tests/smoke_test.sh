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
run_case grafted_comb --strand-topology grafted --backbone-length 64 \
    --side-chain-length 8 --graft-spacing 12 \
    --graft-functional-fraction 40 --strand-count 4 \
    --crosslinker-length 8 --density 0.02 --output data.grafted_comb
run_case grafted_bottlebrush --strand-topology grafted --backbone-length 24 \
    --side-chain-length 5 --graft-spacing 0 \
    --graft-functional-fraction 25 --strand-count 2 \
    --crosslinker-length 8 --density 0.01 \
    --output data.grafted_bottlebrush

base_dir="$test_root/base/PDMS_elastomer"
full_dir="$test_root/full/PDMS_elastomer_filler_N8_5wt_film_H40"
ratio_dir="$test_root/ratio/PDMS_elastomer"
config_dir="$test_root/config/from_config"
ring_bifunctional_dir="$test_root/ring_bifunctional/ring_bifunctional"
ring_tetrafunctional_dir="$test_root/ring_tetrafunctional/ring_tetrafunctional"
star_3arm_dir="$test_root/star_3arm/star_3arm"
star_4arm_dir="$test_root/star_4arm/star_4arm"
star_6arm_dir="$test_root/star_6arm/star_6arm"
star_8arm_dir="$test_root/star_8arm/star_8arm"
grafted_comb_dir="$test_root/grafted_comb/grafted_comb"
grafted_bottlebrush_dir="$test_root/grafted_bottlebrush/grafted_bottlebrush"

test -f "$base_dir/data.PDMS_elastomer"
test -f "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_H40"
test -f "$ratio_dir/data.PDMS_elastomer"
test -f "$config_dir/data.from_config"
test -f "$ring_bifunctional_dir/data.ring_bifunctional"
test -f "$ring_tetrafunctional_dir/data.ring_tetrafunctional"
test -f "$star_3arm_dir/data.star_3arm"
test -f "$star_4arm_dir/data.star_4arm"
test -f "$star_6arm_dir/data.star_6arm"
test -f "$star_8arm_dir/data.star_8arm"
test -f "$grafted_comb_dir/data.grafted_comb"
test -f "$grafted_bottlebrush_dir/data.grafted_bottlebrush"

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
check_data "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_H40"
check_data "$ratio_dir/data.PDMS_elastomer"
check_data "$config_dir/data.from_config"
check_data "$ring_bifunctional_dir/data.ring_bifunctional"
check_data "$ring_tetrafunctional_dir/data.ring_tetrafunctional"
check_data "$star_3arm_dir/data.star_3arm"
check_data "$star_4arm_dir/data.star_4arm"
check_data "$star_6arm_dir/data.star_6arm"
check_data "$star_8arm_dir/data.star_8arm"
check_data "$grafted_comb_dir/data.grafted_comb"
check_data "$grafted_bottlebrush_dir/data.grafted_bottlebrush"

grep -q 'prob 0.10000000' "$base_dir/in.PDMS_elastomer"
grep -q '^comm_modify     cutoff 25$' "$base_dir/in.PDMS_elastomer"
grep -q '^velocity        all create 800.0 5489 mom yes rot yes dist gaussian$' \
    "$base_dir/in.PDMS_elastomer"
test "$(grep -c 'wall/lj126' "$full_dir/in.PDMS_elastomer_filler_N8_5wt_film_H40")" -eq 4
grep -q '^boundary        p p f$' \
    "$full_dir/in.PDMS_elastomer_filler_N8_5wt_film_H40"
grep -q '"film_thickness_source": "requested nominal wall-free thickness at 300 K"' \
    "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"film_thickness_definition": "Lz - 2\*cold_lj.cutoff"' \
    "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"density_definition": "mass / nominal wall-free material volume"' \
    "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
awk '
    FNR == NR {
        if ($1 == "\"film_thickness_angstrom\":") thickness = $2 + 0
        if ($1 == "\"film_wall_cutoff_per_side_angstrom\":") cutoff = $2 + 0
        if ($1 == "\"total_mass_g_per_mol_equivalent\":") mass = $2 + 0
        if ($1 == "\"density_g_per_cm3\":") requested_density = $2 + 0
        next
    }
    /xlo xhi$/ { lx = $2 - $1 }
    /ylo yhi$/ { ly = $2 - $1 }
    /zlo zhi$/ { lz = $2 - $1 }
    END {
        expected_lz = thickness + 2.0 * cutoff
        material_density = mass / (0.602214076 * lx * ly * thickness)
        if (thickness != 40.0 || cutoff <= 0.0 ||
            (lz - expected_lz)^2 > 1.0e-10 ||
            (material_density - requested_density)^2 > 2.5e-9) exit 1
    }
' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info" \
  "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_H40"
grep -q '"format": "pdms-elastomer-model-info"' "$base_dir/PDMS_elastomer.info"
grep -q '"format_version": 3' "$base_dir/PDMS_elastomer.info"
grep -q '"strands": {"component": 1, "N": 16, "M": 12' "$base_dir/PDMS_elastomer.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 6' "$base_dir/PDMS_elastomer.info"
grep -q '"moderators": {"component": 3, "N": 5, "M": 2' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"filler": {"component": 4, "N": 8, "M": 3' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"molecule_id_start": 19, "molecule_id_end": 20' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"molecule_id_start": 21, "molecule_id_end": 23' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"moderators_included": false' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"extra_moderator_functional_groups": 8' "$full_dir/PDMS_elastomer_filler_N8_5wt_film_H40.info"
grep -q '"requested": "2:1"' "$ratio_dir/PDMS_elastomer.info"
grep -q '"crosslinkers": {"component": 2, "N": 16, "M": 3' "$ratio_dir/PDMS_elastomer.info"
grep -q '"expected_frames": 1001' "$base_dir/PDMS_elastomer.info"
grep -q '"initial_velocity_seed": 5489' "$base_dir/PDMS_elastomer.info"
grep -q '"initial_velocity_temperature_K": 800.0' \
    "$base_dir/PDMS_elastomer.info"
grep -q '"hard_maximum_beads": null' "$base_dir/PDMS_elastomer.info"
grep -q '"recommended_maximum_beads": 150000' "$base_dir/PDMS_elastomer.info"
grep -q '"exceeds_recommended_maximum": false' "$base_dir/PDMS_elastomer.info"
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

check_grafted() {
    local directory=$1
    local case_name=$2
    local label=$3
    local backbone=$4
    local side_length=$5
    local spacing=$6
    local side_count=$7
    local functionality=$8
    local strand_beads=$9
    local strand_count=${10}
    local crosslinker_count=${11}
    local data="$directory/data.$case_name"
    local info="$directory/$case_name.info"

    grep -q '"topology": "grafted"' "$info"
    grep -q "\"grafted_label\": \"$label\"" "$info"
    grep -q "\"grafted_backbone_length\": $backbone" "$info"
    grep -q "\"grafted_side_chain_length\": $side_length" "$info"
    grep -q "\"graft_spacing\": $spacing" "$info"
    grep -q "\"side_chain_count\": $side_count" "$info"
    grep -q "\"functional_side_chain_count\": $functionality" "$info"
    grep -q '"reactive_distribution": "selected_side_chain_ends"' "$info"
    grep -q "\"strands\": {\"component\": 1, \"N\": $strand_beads, \"M\": $strand_count" "$info"
    grep -q "\"crosslinkers\": {\"component\": 2, \"N\": 8, \"M\": $crosslinker_count" "$info"

    awk -v molecules="$strand_count" -v beads="$strand_beads" \
        -v backbone="$backbone" -v side_length="$side_length" \
        -v expected="$functionality" '
        /^Atoms # full$/ { in_atoms = 1; next }
        in_atoms && /^Bonds$/ { in_atoms = 0 }
        in_atoms && NF == 7 && $2 <= molecules {
            local = (($1 - 1) % beads) + 1
            if (local <= backbone && $3 != 1) bad_backbone = 1
            if ($3 == 2) {
                ++reactive[$2]
                if (local <= backbone || (local - backbone) % side_length != 0)
                    bad_reactive = 1
            }
        }
        END {
            for (molecule = 1; molecule <= molecules; ++molecule)
                if (reactive[molecule] != expected) exit 1
            if (bad_backbone || bad_reactive) exit 1
        }
    ' "$data"

    awk -v backbone="$backbone" -v expected="$side_count" '
        /^Bonds$/ { in_bonds = 1; next }
        in_bonds && /^[A-Za-z]/ { in_bonds = 0 }
        in_bonds && NF == 4 && (($3 <= backbone && $4 > backbone) ||
            ($4 <= backbone && $3 > backbone)) { ++attachments }
        END { exit attachments != expected }
    ' "$data"
}

check_grafted "$grafted_comb_dir" grafted_comb comb_like 64 8 12 5 2 104 4 2
check_grafted "$grafted_bottlebrush_dir" grafted_bottlebrush \
    dense_bottlebrush 24 5 0 24 6 144 2 3

limit_dir="$test_root/limit"
mkdir -p "$limit_dir"
(
    cd "$limit_dir"
    "$generator" --strand-count 1100 >/dev/null 2>limit.err
)
test -f "$limit_dir/PDMS_elastomer/data.PDMS_elastomer"
grep -q 'warning: total beads exceed the recommended 150000-bead' \
    "$limit_dir/limit.err"
grep -q '"total_beads": 158400' \
    "$limit_dir/PDMS_elastomer/PDMS_elastomer.info"
grep -q '"hard_maximum_beads": null' \
    "$limit_dir/PDMS_elastomer/PDMS_elastomer.info"
grep -q '"exceeds_recommended_maximum": true' \
    "$limit_dir/PDMS_elastomer/PDMS_elastomer.info"

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

invalid_grafted_dir="$test_root/invalid_grafted"
mkdir -p "$invalid_grafted_dir"
for option_value in '--graft-spacing -1' \
    '--graft-functional-fraction 0' \
    '--graft-functional-fraction 101'; do
    read -r option value <<< "$option_value"
    if (
        cd "$invalid_grafted_dir"
        "$generator" --strand-topology grafted "$option" "$value" \
            >/dev/null 2>grafted.err
    ); then
        echo "invalid grafted option unexpectedly succeeded: $option_value" >&2
        exit 1
    fi
done

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
' "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_H40"

topology_analyzer="$repo_root/bin/topology_analyzer"
basic_analyzer="$repo_root/bin/basic_network_analyzer"
profile_analyzer="$repo_root/bin/network_profile_analyzer"

reacted_copy() {
    local source=$1
    local destination=$2
    local pairs=$3
    local additions
    additions=$(awk -F';' '{print NF}' <<< "$pairs")
    awk -v pairs="$pairs" -v additions="$additions" '
        /^[0-9]+ bonds$/ && NF == 2 {
            original_bonds = $1
            print $1 + additions, $2
            next
        }
        /^Bonds$/ { in_bonds = 1; print; next }
        /^Angles$/ && in_bonds {
            split(pairs, entries, ";")
            for (i = 1; i <= additions; ++i) {
                split(entries[i], atoms, ",")
                print original_bonds + i, 2, atoms[1], atoms[2]
            }
            print ""
            in_bonds = 0
        }
        { print }
    ' "$source" > "$destination"
}

mapfile -t ring_xlink_atoms < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $3 == 3 && !seen[$2]++ { print $1 }
' "$ring_bifunctional_dir/data.ring_bifunctional" | head -2)
ring_final="$ring_bifunctional_dir/data.ring_bifunctional.npt_eq"
reacted_copy "$ring_bifunctional_dir/data.ring_bifunctional" "$ring_final" \
    "1,${ring_xlink_atoms[0]};9,${ring_xlink_atoms[1]}"

ring_analysis="$test_root/ring_analysis"
"$topology_analyzer" "$ring_final" \
    "$ring_bifunctional_dir/ring_bifunctional.info" \
    --output-dir "$ring_analysis" --skip-self-paths >/dev/null
grep -q $'status active: 2' \
    "$ring_analysis/topology_report.ring_bifunctional.txt"
grep -q $'fully_reacted_ring' \
    "$ring_analysis/parent_molecules.ring_bifunctional.tsv"
test "$(head -1 "$ring_analysis/config.ring_bifunctional.Z1")" -eq 2
test "$(sed -n '3p' "$ring_analysis/config.ring_bifunctional.Z1")" = "7 7"
test "$(wc -l < "$ring_analysis/config.ring_bifunctional.Z1.map.tsv")" -eq 3
awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; ++column) header[$column] = column
        next
    }
    {
        if ($header["effective_contour_beads"] != 9 ||
            $header["z1_beads"] != 7) exit 1
        split($header["excluded_atom_ids"], excluded, ",")
        for (i in excluded) excluded_atoms[excluded[i]] = 1
        count = split($header["z1_atom_ids"], atoms, ",")
        if (count != 7) exit 1
        for (i = 1; i <= count; ++i) {
            if (atoms[i] == 1 || atoms[i] == 9 || seen[atoms[i]]++) exit 1
        }
        ++rows
    }
    END { exit rows != 2 || !excluded_atoms[1] || !excluded_atoms[9] }
' "$ring_analysis/config.ring_bifunctional.Z1.map.tsv"
read -r z1_lx z1_ly z1_lz < <(sed -n '2p' \
    "$ring_analysis/config.ring_bifunctional.Z1")
read -r data_lx data_ly data_lz < <(awk '
    /xlo xhi$/ { lx = $2 - $1 }
    /ylo yhi$/ { ly = $2 - $1 }
    /zlo zhi$/ { lz = $2 - $1 }
    END { printf "%.15g %.15g %.15g\n", lx, ly, lz }
' "$ring_final")
awk -v z1_lx="$z1_lx" -v z1_ly="$z1_ly" -v z1_lz="$z1_lz" \
    -v data_lx="$data_lx" -v data_ly="$data_ly" -v data_lz="$data_lz" '
    function abs(value) { return value < 0 ? -value : value }
    BEGIN {
        if (abs(z1_lx - data_lx) > 1e-9 ||
            abs(z1_ly - data_ly) > 1e-9 ||
            abs(z1_lz - data_lz) > 1e-9) exit 1
    }
'
grep -q $'\t1$' "$ring_analysis/config.ring_bifunctional.Z1.map.tsv"

mapfile -t ring_same_xlink_atoms < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $2 == 13 && $3 == 3 { print $1 }
' "$ring_bifunctional_dir/data.ring_bifunctional" | head -2)

ring_dangling_final="$ring_bifunctional_dir/data.ring_bifunctional.dangling.npt_eq"
reacted_copy "$ring_bifunctional_dir/data.ring_bifunctional" \
    "$ring_dangling_final" "1,${ring_same_xlink_atoms[0]}"
ring_dangling_analysis="$test_root/ring_dangling_analysis"
"$topology_analyzer" "$ring_dangling_final" \
    "$ring_bifunctional_dir/ring_bifunctional.info" \
    --output-dir "$ring_dangling_analysis" --skip-self-paths >/dev/null
test "$(head -1 "$ring_dangling_analysis/config.ring_bifunctional.Z1")" -eq 1
test "$(sed -n '3p' "$ring_dangling_analysis/config.ring_bifunctional.Z1")" = "15"
awk -F '\t' 'NR == 2 { exit !($5 == "dangling_loop" && $10 == 17 && $11 == 15) }' \
    "$ring_dangling_analysis/config.ring_bifunctional.Z1.map.tsv"

ring_self_final="$ring_bifunctional_dir/data.ring_bifunctional.self_loop.npt_eq"
reacted_copy "$ring_bifunctional_dir/data.ring_bifunctional" \
    "$ring_self_final" \
    "1,${ring_same_xlink_atoms[0]};9,${ring_same_xlink_atoms[1]}"
ring_self_analysis="$test_root/ring_self_analysis"
"$topology_analyzer" "$ring_self_final" \
    "$ring_bifunctional_dir/ring_bifunctional.info" \
    --output-dir "$ring_self_analysis" --skip-self-paths >/dev/null
test "$(head -1 "$ring_self_analysis/config.ring_bifunctional.Z1")" -eq 2
test "$(sed -n '3p' "$ring_self_analysis/config.ring_bifunctional.Z1")" = "7 7"
test "$(awk -F '\t' 'NR > 1 && $5 == "self_loop" { ++count } END { print count + 0 }' \
    "$ring_self_analysis/config.ring_bifunctional.Z1.map.tsv")" -eq 2

ring_final_velocities="$ring_bifunctional_dir/data.ring_bifunctional.npt_eq.velocities"
awk '
    /^[0-9]+ atoms$/ && NF == 2 { atom_count = $1 }
    { print }
    END {
        print ""
        print "Velocities"
        print ""
        for (id = 1; id <= atom_count; ++id) {
            vx = id % 2 ? 0.01 : -0.01
            print id, vx, 0.0, 0.0
        }
    }
' "$ring_final" > "$ring_final_velocities"

basic_analysis="$test_root/basic_analysis"
"$basic_analyzer" "$ring_final_velocities" \
    "$ring_bifunctional_dir/ring_bifunctional.info" \
    --histogram-bins 8 --output-dir "$basic_analysis" >/dev/null
for output in \
    "basic_network_report.ring_bifunctional.txt" \
    "network_statistics.ring_bifunctional.tsv" \
    "strand_properties.ring_bifunctional.tsv" \
    "strand_statistics.ring_bifunctional.tsv" \
    "strand_histograms.ring_bifunctional.tsv" \
    "junction_properties.ring_bifunctional.tsv"; do
    test -s "$basic_analysis/$output"
done
test "$(wc -l < "$basic_analysis/strand_properties.ring_bifunctional.tsv")" -eq 3
awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; ++column) header[$column] = column
        next
    }
    {
        if ($header["status"] != "active" || $header["Ree_A"] <= 0 ||
            $header["Rg_A"] <= 0 || $header["Lpp_A"] != "nan" ||
            $header["Lpp_source"] != "not_available_without_primitive_path_analysis")
            exit 1
        ++rows
    }
    END { exit rows != 2 }
' "$basic_analysis/strand_properties.ring_bifunctional.tsv"
awk -F '\t' '
    $1 == "active" && $2 == "Ree_A" && $3 == 2 && $4 > 0 { ree = 1 }
    $1 == "active" && $2 == "Rg_A" && $3 == 2 && $4 > 0 { rg = 1 }
    $1 == "active" && $2 == "Lpp_A" && $3 == 0 && $4 == "nan" { lpp = 1 }
    END { exit !(ree && rg && lpp) }
' "$basic_analysis/strand_statistics.ring_bifunctional.tsv"
awk -F '\t' '
    $1 == "density" && $2 > 0 { density = 1 }
    $1 == "velocity_temperature" && $2 > 0 { temperature = 1 }
    $1 == "active_strands" && $2 == 2 { active = 1 }
    $1 == "Lpp_available" && $2 == 0 { lpp = 1 }
    END { exit !(density && temperature && active && lpp) }
' "$basic_analysis/network_statistics.ring_bifunctional.tsv"
grep -q 'Lpp: unavailable (NaN)' \
    "$basic_analysis/basic_network_report.ring_bifunctional.txt"

mapfile -t star_xlink_atoms < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $3 == 3 && !seen[$2]++ { print $1 }
' "$star_4arm_dir/data.star_4arm" | head -4)
star_final="$star_4arm_dir/data.star_4arm.npt_eq"
reacted_copy "$star_4arm_dir/data.star_4arm" "$star_final" \
    "5,${star_xlink_atoms[0]};9,${star_xlink_atoms[1]};13,${star_xlink_atoms[2]};17,${star_xlink_atoms[3]}"
star_analysis="$test_root/star_analysis"
"$topology_analyzer" "$star_final" "$star_4arm_dir/star_4arm.info" \
    --output-dir "$star_analysis" --skip-self-paths >/dev/null
grep -q 'status active: 4' "$star_analysis/topology_report.star_4arm.txt"
grep -q $'star_center' "$star_analysis/network_nodes.star_4arm.tsv"
test "$(head -1 "$star_analysis/config.star_4arm.Z1")" -eq 4
test "$(sed -n '3p' "$star_analysis/config.star_4arm.Z1")" = "4 4 4 4"
awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; ++column) header[$column] = column
        next
    }
    {
        if ($header["effective_contour_beads"] != 5 ||
            $header["z1_beads"] != 4 ||
            $header["excluded_atom_ids"] != 1) exit 1
        count = split($header["z1_atom_ids"], atoms, ",")
        if (count != 4) exit 1
        for (i = 1; i <= count; ++i)
            if (atoms[i] == 1 || seen[atoms[i]]++) exit 1
        ++rows
    }
    END { exit rows != 4 }
' "$star_analysis/config.star_4arm.Z1.map.tsv"

star_partial_final="$star_4arm_dir/data.star_4arm.partial.npt_eq"
reacted_copy "$star_4arm_dir/data.star_4arm" "$star_partial_final" \
    "5,${star_xlink_atoms[0]}"
star_partial_analysis="$test_root/star_partial_analysis"
"$topology_analyzer" "$star_partial_final" "$star_4arm_dir/star_4arm.info" \
    --output-dir "$star_partial_analysis" --skip-self-paths >/dev/null
test "$(head -1 "$star_partial_analysis/config.star_4arm.Z1")" -eq 4
test "$(awk -F '\t' 'NR > 1 && $5 == "active" { ++count } END { print count + 0 }' \
    "$star_partial_analysis/config.star_4arm.Z1.map.tsv")" -eq 1
test "$(awk -F '\t' 'NR > 1 && $5 == "dangling" { ++count } END { print count + 0 }' \
    "$star_partial_analysis/config.star_4arm.Z1.map.tsv")" -eq 3

mapfile -t graft_sites < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $2 == 1 && $3 == 2 { print $1 }
' "$grafted_comb_dir/data.grafted_comb")
mapfile -t graft_xlink_atoms < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $3 == 3 && !seen[$2]++ { print $1 }
' "$grafted_comb_dir/data.grafted_comb" | head -2)
graft_final="$grafted_comb_dir/data.grafted_comb.npt_eq"
reacted_copy "$grafted_comb_dir/data.grafted_comb" "$graft_final" \
    "${graft_sites[0]},${graft_xlink_atoms[0]};${graft_sites[1]},${graft_xlink_atoms[1]}"
graft_analysis="$test_root/graft_analysis"
"$topology_analyzer" "$graft_final" "$grafted_comb_dir/grafted_comb.info" \
    --output-dir "$graft_analysis" --skip-self-paths >/dev/null
grep -q 'status active: 1' "$graft_analysis/topology_report.grafted_comb.txt"
test "$(head -1 "$graft_analysis/config.grafted_comb.Z1")" -eq 2
test "$(sed -n '3p' "$graft_analysis/config.grafted_comb.Z1")" = "34 8"
awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; ++column) header[$column] = column
        next
    }
    {
        count = split($header["z1_atom_ids"], atoms, ",")
        if (count != $header["z1_beads"]) exit 1
        for (i = 1; i <= count; ++i)
            if (seen[atoms[i]]++) exit 1
        if (NR == 2) {
            if ($header["effective_strand_id"] == 0 ||
                $header["status"] != "active" ||
                $header["effective_contour_beads"] != 43 ||
                $header["z1_beads"] != 34 ||
                $header["reacted_site_index"] != 1 ||
                $header["first_atom"] != 80 ||
                $header["second_atom"] != 39 ||
                $header["next_graft_atom"] != 40) exit 1
        } else if (NR == 3) {
            if ($header["effective_strand_id"] != 0 ||
                $header["status"] != "dangling" ||
                $header["effective_contour_beads"] != 9 ||
                $header["z1_beads"] != 8 ||
                $header["reacted_site_index"] != 2 ||
                $header["first_atom"] != 96 ||
                $header["second_atom"] != 89 ||
                $header["next_graft_atom"] != 0) exit 1
        }
        ++rows
    }
    END { exit rows != 2 }
' "$graft_analysis/config.grafted_comb.Z1.map.tsv"

graft_single_final="$grafted_comb_dir/data.grafted_comb.single.npt_eq"
reacted_copy "$grafted_comb_dir/data.grafted_comb" "$graft_single_final" \
    "${graft_sites[0]},${graft_xlink_atoms[0]}"
graft_single_analysis="$test_root/graft_single_analysis"
"$topology_analyzer" "$graft_single_final" \
    "$grafted_comb_dir/grafted_comb.info" \
    --output-dir "$graft_single_analysis" --skip-self-paths >/dev/null
test "$(head -1 "$graft_single_analysis/config.grafted_comb.Z1")" -eq 1
test "$(sed -n '3p' "$graft_single_analysis/config.grafted_comb.Z1")" = "8"
awk -F '\t' '
    NR == 1 { for (column = 1; column <= NF; ++column) header[$column] = column; next }
    NR == 2 {
        exit !($header["effective_strand_id"] == 0 &&
            $header["status"] == "dangling" &&
            $header["reacted_site_index"] == 1 &&
            $header["first_atom"] == 80 && $header["second_atom"] == 73 &&
            $header["z1_beads"] == 8)
    }
' "$graft_single_analysis/config.grafted_comb.Z1.map.tsv"

mapfile -t graft_same_xlink_atoms < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $2 == 5 && $3 == 3 { print $1 }
' "$grafted_comb_dir/data.grafted_comb" | head -2)
graft_self_final="$grafted_comb_dir/data.grafted_comb.self_loop.npt_eq"
reacted_copy "$grafted_comb_dir/data.grafted_comb" "$graft_self_final" \
    "${graft_sites[0]},${graft_same_xlink_atoms[0]};${graft_sites[1]},${graft_same_xlink_atoms[1]}"
graft_self_analysis="$test_root/graft_self_analysis"
"$topology_analyzer" "$graft_self_final" \
    "$grafted_comb_dir/grafted_comb.info" \
    --output-dir "$graft_self_analysis" --skip-self-paths >/dev/null
test "$(head -1 "$graft_self_analysis/config.grafted_comb.Z1")" -eq 2
test "$(awk -F '\t' 'NR > 1 && $5 == "self_loop" { ++count } END { print count + 0 }' \
    "$graft_self_analysis/config.grafted_comb.Z1.map.tsv")" -eq 1
test "$(awk -F '\t' 'NR > 1 && $5 == "dangling" { ++count } END { print count + 0 }' \
    "$graft_self_analysis/config.grafted_comb.Z1.map.tsv")" -eq 1

mapfile -t bottle_graft_sites < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $2 == 1 && $3 == 2 { print $1 }
' "$grafted_bottlebrush_dir/data.grafted_bottlebrush")
mapfile -t bottle_xlink_atoms < <(awk '
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { exit }
    in_atoms && NF >= 7 && $3 == 3 && !seen[$2]++ { print $1 }
' "$grafted_bottlebrush_dir/data.grafted_bottlebrush" | head -2)
bottle_final="$grafted_bottlebrush_dir/data.grafted_bottlebrush.npt_eq"
reacted_copy "$grafted_bottlebrush_dir/data.grafted_bottlebrush" \
    "$bottle_final" \
    "${bottle_graft_sites[0]},${bottle_xlink_atoms[0]};${bottle_graft_sites[1]},${bottle_xlink_atoms[1]}"
bottle_analysis="$test_root/bottle_analysis"
"$topology_analyzer" "$bottle_final" \
    "$grafted_bottlebrush_dir/grafted_bottlebrush.info" \
    --output-dir "$bottle_analysis" --skip-self-paths >/dev/null
test "$(head -1 "$bottle_analysis/config.grafted_bottlebrush.Z1")" -eq 2
test "$(sed -n '3p' "$bottle_analysis/config.grafted_bottlebrush.Z1")" = "9 5"
awk -F '\t' '
    NR == 1 { for (column = 1; column <= NF; ++column) header[$column] = column; next }
    {
        count = split($header["z1_atom_ids"], atoms, ",")
        for (i = 1; i <= count; ++i) if (seen[atoms[i]]++) exit 1
        ++rows
    }
    END { exit rows != 2 }
' "$bottle_analysis/config.grafted_bottlebrush.Z1.map.tsv"

trajectory="$test_root/ring_two_frames.lammpstrj"
awk '
    /^[0-9]+ atoms$/ && NF == 2 { atom_count = $1 }
    /xlo xhi$/ { xlo = $1; xhi = $2 }
    /ylo yhi$/ { ylo = $1; yhi = $2 }
    /zlo zhi$/ { zlo = $1; zhi = $2 }
    /^Atoms # full$/ { in_atoms = 1; next }
    in_atoms && /^Bonds$/ { in_atoms = 0 }
    in_atoms && NF >= 7 {
        id[++count] = $1; molecule[count] = $2; type[count] = $3
        x[count] = $5; y[count] = $6; z[count] = $7
    }
    END {
        for (frame = 0; frame < 2; ++frame) {
            print "ITEM: TIMESTEP"
            print frame * 1000
            print "ITEM: NUMBER OF ATOMS"
            print atom_count
            print "ITEM: BOX BOUNDS pp pp pp"
            print xlo, xhi; print ylo, yhi; print zlo, zhi
            print "ITEM: ATOMS id mol type x y z ix iy iz"
            for (i = 1; i <= count; ++i) {
                displacement = frame && molecule[i] <= 12 ? 1.0 : 0.0
                print id[i], molecule[i], type[i], x[i] + displacement,
                      y[i], z[i], 0, 0, 0
            }
        }
    }
' "$ring_final" > "$trajectory"

profile_analysis="$test_root/profile_analysis"
"$profile_analyzer" "$ring_final" \
    "$ring_bifunctional_dir/ring_bifunctional.info" \
    --trajectory "$trajectory" --bin-width 20 \
    --output-dir "$profile_analysis" >/dev/null
test -s "$profile_analysis/network_z_profile.ring_bifunctional.tsv"
test -s "$profile_analysis/layer_dynamics.ring_bifunctional.tsv"
awk 'NR > 1 && $1 == 1 && $11 > 0 { found = 1 } END { exit !found }' \
    "$profile_analysis/layer_dynamics.ring_bifunctional.tsv"

echo "PDMS elastomer smoke tests passed"
