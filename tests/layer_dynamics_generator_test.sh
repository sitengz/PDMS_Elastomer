#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
model_generator="$repo_root/bin/pdms_elastomer_generator"
dynamics_generator="$repo_root/bin/layer_dynamics_generator"
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
    awk '
        { print }
        /^Atoms # full$/ {
            saw_atoms = 1
        }
        END {
            if (!saw_atoms) exit 1
            print ""
            print "Velocities"
            print ""
            for (i = 1; i <= 128; ++i) print i, 0.0, 0.0, 0.0
        }
    ' "$test_root/$name/data.$name.npt_eq" > "$test_root/$name/data.$name.npt_eq.with_velocities"
    mv "$test_root/$name/data.$name.npt_eq.with_velocities" \
       "$test_root/$name/data.$name.npt_eq"
}

generate_source dynamics_bulk
generate_source dynamics_film --thickness 40

"$dynamics_generator" \
    "$test_root/dynamics_bulk/data.dynamics_bulk.npt_eq" \
    "$test_root/dynamics_bulk/dynamics_bulk.info" >/dev/null
"$dynamics_generator" \
    "$test_root/dynamics_film/data.dynamics_film.npt_eq" \
    "$test_root/dynamics_film/dynamics_film.info" >/dev/null

bulk_input="$test_root/dynamics_bulk/layer_dynamics/in.layer_dynamics.dynamics_bulk"
film_input="$test_root/dynamics_film/layer_dynamics/in.layer_dynamics.dynamics_film"

for input in "$bulk_input" "$film_input"; do
    test -f "$input"
    grep -q '^comm_modify     cutoff 25$' "$input"
    grep -q '^pair_style      lj/gromacs 12 15$' "$input"
    test "$(grep -c '^reset_timestep  0$' "$input")" -eq 2
    grep -q '^compute         global_dw all msd com yes$' "$input"
    grep -q '^compute         global_msd all msd com yes$' "$input"
    grep -q '^dump            dw all custom 20 dump.debye_waller.' "$input"
    grep -q '^dump_modify     dw first yes sort id$' "$input"
    grep -q '^run             20000$' "$input"
    grep -q '^undump          dw$' "$input"
    grep -q '^dump            msd all custom 1000 dump.layer_dynamics.' "$input"
    grep -q '^dump_modify     msd first yes sort id$' "$input"
    grep -q '^dump_modify     msd every 5000 first no$' "$input"
    test "$(grep -c '^run             1000000$' "$input")" -eq 2
    grep -q '^run             4000000$' "$input"
    grep -q '^fix             dw_integrate all nvt temp 300.000000000 300.000000000 50.0$' "$input"
    grep -q '^fix             msd_integrate all nvt temp 300.000000000 300.000000000 50.0$' "$input"
    if grep -Eq '^(velocity|minimize|fix +xlink|fix +deform_box)' "$input"; then
        echo "Layer-dynamics input unexpectedly changes the source state or topology" >&2
        exit 1
    fi
done

grep -q '^boundary        p p p$' "$bulk_input"
test "$(grep -c 'wall/lj126' "$bulk_input")" -eq 0
grep -q '^boundary        p p f$' "$film_input"
test "$(grep -c 'wall/lj126' "$film_input")" -eq 0
grep -q '^change_box      all z delta -50.000000000 50.000000000 units box$' \
    "$film_input"
grep -q '^write_data      data.dynamics_film.free_surface_eq nocoeff$' "$film_input"

film_info="$test_root/dynamics_film/layer_dynamics/layer_dynamics.dynamics_film.info"
grep -q '"surface_equilibration_steps": 1000000' "$film_info"
grep -q '"debye_waller_steps": 20000' "$film_info"
grep -q '"debye_waller_duration_ps": 100.0000000000' "$film_info"
grep -q '"debye_waller_dump_every_ps": 0.1000000000' "$film_info"
grep -q '"expected_debye_waller_frames": 1001' "$film_info"
grep -q '"production_steps": 5000000' "$film_info"
grep -q '"production_duration_ns": 25.0000000000' "$film_info"
grep -q '"expected_msd_trajectory_frames": 1801' "$film_info"
grep -q '"source_velocities_retained": true' "$film_info"
grep -q '"each_trajectory_first_frame_is_origin": true' "$film_info"
grep -q '"source_walls_recreated": false' "$film_info"
grep -q '"free_surfaces": true' "$film_info"
grep -q '"vacuum_padding_per_side_angstrom": 50.0000000000' "$film_info"

test -x "$test_root/dynamics_film/layer_dynamics/submit.layer_dynamics.dynamics_film.sh"

override_output="$test_root/wrapper_output"
NPT_EQ_FILE="$test_root/dynamics_film/data.dynamics_film.npt_eq" \
INFO_FILE="$test_root/dynamics_film/dynamics_film.info" \
LAYER_DYNAMICS_OUTPUT_DIR="$override_output" \
    bash "$repo_root/simulations/01_linear_reference/film_2Ree/run_layer_dynamics.sh" \
        --production-steps 2000000 >/dev/null
test -f "$override_output/in.layer_dynamics.dynamics_film"
grep -q '^run             1000000$' "$override_output/in.layer_dynamics.dynamics_film"

echo "Layer-dynamics generator tests passed"
