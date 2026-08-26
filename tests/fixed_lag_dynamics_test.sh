#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
model_generator="$repo_root/bin/pdms_elastomer_generator"
profile_analyzer="$repo_root/bin/network_profile_analyzer"
fixed_analyzer="$repo_root/bin/fixed_lag_dynamics_analyzer"
test_root=$(mktemp -d)
trap 'rm -rf -- "$test_root"' EXIT

(
    cd "$test_root"
    "$model_generator" --strand-length 12 --strand-count 8 \
        --crosslinker-length 8 --output data.fixed_lag >/dev/null
)
case_dir="$test_root/fixed_lag"
data_file="$case_dir/data.fixed_lag.npt_eq"
info_file="$case_dir/fixed_lag.info"
cp "$case_dir/data.fixed_lag" "$data_file"

make_dump() {
    local output=$1
    local step_increment=$2
    awk -v step_increment="$step_increment" '
        / atoms$/ { atom_count = $1 }
        /xlo xhi$/ { xlo = $1; xhi = $2 }
        /ylo yhi$/ { ylo = $1; yhi = $2 }
        /zlo zhi$/ { zlo = $1; zhi = $2 }
        /^Atoms # full$/ { in_atoms = 1; next }
        /^Bonds$/ { in_atoms = 0 }
        in_atoms && NF == 7 && $1 ~ /^[0-9]+$/ {
            id[$1] = $1
            molecule[$1] = $2
            x[$1] = $5
            y[$1] = $6
            z[$1] = $7
        }
        END {
            if (atom_count <= 0 || length(id) != atom_count) exit 1
            for (frame = 0; frame < 5; ++frame) {
                print "ITEM: TIMESTEP"
                print frame * step_increment
                print "ITEM: NUMBER OF ATOMS"
                print atom_count
                print "ITEM: BOX BOUNDS pp pp pp"
                print xlo, xhi
                print ylo, yhi
                print zlo, zhi
                print "ITEM: ATOMS id xu yu zu"
                for (atom = 1; atom <= atom_count; ++atom) {
                    displacement = molecule[atom] <= 8 ? 0.2 * frame : 0.0
                    print atom, x[atom] + displacement, y[atom], z[atom]
                }
            }
        }
    ' "$data_file" > "$output"
}

dw_dump="$case_dir/dump.debye_waller.fixed_lag.lammpstrj"
msd_dump="$case_dir/dump.layer_dynamics.fixed_lag.lammpstrj"
make_dump "$dw_dump" 1000
make_dump "$msd_dump" 1000000

profile_output="$case_dir/profile_output"
"$profile_analyzer" "$data_file" "$info_file" \
    --dw-trajectory "$dw_dump" --time-averaged-dw \
    --dw-time-ps 10 --dw-time-origin-count 3 \
    --dw-time-origin-stride 1 --output-dir "$profile_output" >/dev/null

test -f "$profile_output/debye_waller_global.fixed_lag.tsv"
test -f "$profile_output/debye_waller_time_averaged.fixed_lag.tsv"
test -f "$profile_output/layer_debye_waller_time_averaged.fixed_lag.tsv"
awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; ++i) {
            if ($i == "time_ps") time_column = i
            if ($i == "time_origins") origins_column = i
            if ($i == "selected") selected_column = i
        }
        next
    }
    $selected_column == 1 {
        if (($time_column - 10)^2 > 1e-12 || $origins_column != 3) exit 1
        found = 1
    }
    END { exit !found }
' "$profile_output/debye_waller_time_averaged.fixed_lag.tsv"

fixed_output="$case_dir/fixed_output"
"$fixed_analyzer" "$data_file" "$info_file" \
    --msd-trajectory "$msd_dump" --dw-trajectory "$dw_dump" \
    --output-dir "$fixed_output" >/dev/null

for output in \
    fixed_lag_msd_origins.fixed_lag.tsv \
    fixed_lag_msd_layers.fixed_lag.tsv \
    fixed_lag_msd_layer_summary.fixed_lag.tsv \
    fixed_lag_msd_layer_pooled_summary.fixed_lag.tsv \
    fixed_lag_u2_origins.fixed_lag.tsv \
    fixed_lag_u2_layers.fixed_lag.tsv \
    fixed_lag_u2_layer_summary.fixed_lag.tsv \
    fixed_lag_u2_layer_pooled_summary.fixed_lag.tsv \
    fixed_lag_summary.fixed_lag.tsv \
    fixed_lag_report.fixed_lag.txt
do
    test -f "$fixed_output/$output"
done

test "$(wc -l < "$fixed_output/fixed_lag_msd_origins.fixed_lag.tsv")" -eq 4
test "$(wc -l < "$fixed_output/fixed_lag_u2_origins.fixed_lag.tsv")" -eq 4
awk -F '\t' '
    NR == 1 { next }
    $1 == "msd_fixed_lag" {
        if ($2 != 10000 || $3 != 3 || ($4 - 10000)^2 > 1e-12) exit 1
        msd = 1
    }
    $1 == "u2_fixed_lag" {
        if ($2 != 10 || $3 != 3 || ($4 - 10)^2 > 1e-12) exit 1
        dw = 1
    }
    END { exit !(msd && dw) }
' "$fixed_output/fixed_lag_summary.fixed_lag.tsv"
awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; ++i) {
            if ($i == "actual_lag_ps") lag = i
            if ($i == "msd_xy_A2") value = i
            if ($i == "D_E_xy_A2_per_ps") diffusion = i
        }
        next
    }
    {
        expected = $value / (4.0 * $lag)
        if (($diffusion - expected)^2 > 1e-24) exit 1
    }
' "$fixed_output/fixed_lag_msd_origins.fixed_lag.tsv"
grep -q 'origins may be correlated' \
    "$fixed_output/fixed_lag_report.fixed_lag.txt"
grep -q 'bead-origin weighted' \
    "$fixed_output/fixed_lag_report.fixed_lag.txt"

raw_layers="$fixed_output/fixed_lag_msd_layers.fixed_lag.tsv"
pooled_layers="$fixed_output/fixed_lag_msd_layer_pooled_summary.fixed_lag.tsv"
awk -F '\t' '
    FNR == NR {
        if (FNR == 1) {
            for (i = 1; i <= NF; ++i) {
                if ($i == "bin") raw_bin = i
                if ($i == "strand_beads") raw_beads = i
                if ($i == "msd_xy_A2") raw_value = i
            }
            next
        }
        weighted[$raw_bin] += $raw_beads * $raw_value
        observations[$raw_bin] += $raw_beads
        next
    }
    FNR == 1 {
        for (i = 1; i <= NF; ++i) {
            if ($i == "bin") pooled_bin = i
            if ($i == "strand_bead_observations") pooled_beads = i
            if ($i == "pooled_mean_xy_A2") pooled_value = i
            if ($i == "pooled_D_E_xy_A2_per_ps") pooled_diffusion = i
        }
        next
    }
    $pooled_beads > 0 {
        expected = weighted[$pooled_bin] / observations[$pooled_bin]
        if (($pooled_value - expected)^2 > 1e-20) exit 1
        if (($pooled_diffusion - expected / 40000.0)^2 > 1e-24) exit 1
        checked = 1
    }
    END { exit !checked }
' "$raw_layers" "$pooled_layers"

declared_atoms=$(awk '$2 == "atoms" { print $1; exit }' "$data_file")
awk -F '\t' -v declared_atoms="$declared_atoms" '
    NR == 1 {
        for (i = 1; i <= NF; ++i)
            if ($i == "mean_all_beads") all_beads = i
        next
    }
    { total += $all_beads }
    END { if ((total - declared_atoms)^2 > 1e-12) exit 1 }
' "$pooled_layers"

echo "Fixed-lag dynamics analyzer tests passed"
