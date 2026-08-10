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
    --crosslinker-length 16 --spacing 3
run_case full --filler-wt 5 --filler-length 8 \
    --strand-length 24 --strand-count 12 --crosslinker-length 16 \
    --moderator-count 2 --spacing 3 \
    --thickness 40
run_case ratio --strand-length 16 --strand-count 12 \
    --crosslinker-length 16 --stoichiometry 2:1 --spacing 3
run_case config --config "$repo_root/tests/basic_model.conf" --strand-count 16

base_dir="$test_root/base/PDMS_elastomer"
full_dir="$test_root/full/PDMS_elastomer_filler_N8_5wt_film_Lz40"
ratio_dir="$test_root/ratio/PDMS_elastomer"
config_dir="$test_root/config/from_config"

test -f "$base_dir/data.PDMS_elastomer"
test -f "$full_dir/data.PDMS_elastomer_filler_N8_5wt_film_Lz40"
test -f "$ratio_dir/data.PDMS_elastomer"
test -f "$config_dir/data.from_config"

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
