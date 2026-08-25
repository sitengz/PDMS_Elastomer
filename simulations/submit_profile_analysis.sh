#!/usr/bin/env bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=nova
#SBATCH --mail-user=siteng@iastate.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

script_path=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")
repo_root=$(cd -- "$(dirname -- "$script_path")/.." && pwd)
analyzer="$repo_root/bin/network_profile_analyzer"

usage() {
    cat <<'EOF'
Usage:
  bash simulations/submit_profile_analysis.sh [options] [SEARCH_ROOT ...]

Discover completed layer-dynamics trajectories below each search root and
submit one network-profile analysis job per case. With no search root, the
repository simulations/ directory is searched.

Options:
  --time-origin-count N   maximum time origins per case (default: 10)
  --time-origin-stride N  origin spacing in analyzed frames (default: 100)
  --frame-stride N        analyze every Nth dump frame (default: 1)
  --bin-width X           target z-bin width in angstrom (default: 5)
  --dry-run               print sbatch commands without submitting
  --help                  show this help

Examples:
  bash simulations/submit_profile_analysis.sh simulations/01_linear_reference

  bash simulations/submit_profile_analysis.sh \
      --time-origin-count 20 --time-origin-stride 50 simulations

Expected case layout:
  CASE_DIR/data.CASE.npt_eq
  CASE_DIR/CASE.info
  CASE_DIR/layer_dynamics/dump.debye_waller.CASE.lammpstrj
  CASE_DIR/layer_dynamics/dump.layer_dynamics.CASE.lammpstrj
EOF
}

require_positive_integer() {
    local value=$1
    local option=$2
    if [[ ! $value =~ ^[1-9][0-9]*$ ]]; then
        echo "$option must be a positive integer: $value" >&2
        exit 2
    fi
}

require_positive_number() {
    local value=$1
    local option=$2
    if [[ ! $value =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] ||
       ! awk -v value="$value" 'BEGIN { exit !(value > 0) }'; then
        echo "$option must be a positive number: $value" >&2
        exit 2
    fi
}

run_worker() {
    if [[ $# -ne 10 ]]; then
        echo "Internal worker argument error" >&2
        exit 2
    fi
    local data_file=$1
    local info_file=$2
    local trajectory_file=$3
    local dw_trajectory_file=$4
    local output_directory=$5
    local origin_count=$6
    local origin_stride=$7
    local frame_stride=$8
    local bin_width=$9
    local case_name=${10}

    if [[ ! -x $analyzer ]]; then
        echo "Missing analyzer executable: $analyzer" >&2
        echo "Run make bin/network_profile_analyzer before submitting." >&2
        exit 2
    fi
    for required_file in "$data_file" "$info_file" "$trajectory_file"; do
        if [[ ! -f $required_file ]]; then
            echo "Missing required analysis input: $required_file" >&2
            exit 2
        fi
    done

    mkdir -p -- "$output_directory"
    local -a command=(
        "$analyzer"
        "$data_file"
        "$info_file"
        --trajectory "$trajectory_file"
        --time-averaged-msd
        --time-origin-count "$origin_count"
        --time-origin-stride "$origin_stride"
        --frame-stride "$frame_stride"
        --bin-width "$bin_width"
        --output-dir "$output_directory"
    )
    if [[ -n $dw_trajectory_file ]]; then
        if [[ ! -f $dw_trajectory_file ]]; then
            echo "Missing Debye-Waller trajectory: $dw_trajectory_file" >&2
            exit 2
        fi
        command+=(--dw-trajectory "$dw_trajectory_file")
    fi

    echo "Case: $case_name"
    echo "Data: $data_file"
    echo "Info: $info_file"
    echo "Layer trajectory: $trajectory_file"
    if [[ -n $dw_trajectory_file ]]; then
        echo "Debye-Waller trajectory: $dw_trajectory_file"
    else
        echo "Debye-Waller trajectory: not found; skipping DW analysis"
    fi
    echo "Output: $output_directory"
    echo "Time origins: at most $origin_count, stride $origin_stride analyzed frames"
    echo "Frame stride: $frame_stride"
    echo "Bin width: $bin_width A"
    echo

    export OMP_NUM_THREADS=1
    srun --ntasks=1 "${command[@]}"
}

if [[ ${1:-} == "--worker" ]]; then
    shift
    run_worker "$@"
    exit 0
fi

time_origin_count=10
time_origin_stride=100
frame_stride=1
bin_width=5
dry_run=0
declare -a search_roots=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --time-origin-count)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            time_origin_count=$2
            shift 2
            ;;
        --time-origin-stride)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            time_origin_stride=$2
            shift 2
            ;;
        --frame-stride)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            frame_stride=$2
            shift 2
            ;;
        --bin-width)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            bin_width=$2
            shift 2
            ;;
        --dry-run)
            dry_run=1
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        --*)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            search_roots+=("$1")
            shift
            ;;
    esac
done

require_positive_integer "$time_origin_count" "--time-origin-count"
if (( time_origin_count < 2 )); then
    echo "--time-origin-count must be at least 2" >&2
    exit 2
fi
require_positive_integer "$time_origin_stride" "--time-origin-stride"
require_positive_integer "$frame_stride" "--frame-stride"
require_positive_number "$bin_width" "--bin-width"

if [[ ${#search_roots[@]} -eq 0 ]]; then
    search_roots=("$repo_root/simulations")
fi

for index in "${!search_roots[@]}"; do
    if [[ ! -d ${search_roots[index]} ]]; then
        echo "Search root is not a directory: ${search_roots[index]}" >&2
        exit 2
    fi
    search_roots[index]=$(cd -- "${search_roots[index]}" && pwd)
done

make -C "$repo_root" bin/network_profile_analyzer

declare -A seen_trajectories=()
submitted=0
skipped=0
for search_root in "${search_roots[@]}"; do
    while IFS= read -r -d '' trajectory_file; do
        if [[ -n ${seen_trajectories[$trajectory_file]+present} ]]; then
            continue
        fi
        seen_trajectories[$trajectory_file]=1

        trajectory_name=$(basename -- "$trajectory_file")
        case_name=${trajectory_name#dump.layer_dynamics.}
        case_name=${case_name%.lammpstrj}
        layer_directory=$(dirname -- "$trajectory_file")
        case_directory=$(dirname -- "$layer_directory")
        data_file="$case_directory/data.$case_name.npt_eq"
        info_file="$case_directory/$case_name.info"
        dw_trajectory_file="$layer_directory/dump.debye_waller.$case_name.lammpstrj"
        output_directory="$case_directory/analysis_$case_name"

        if [[ ! -f $data_file || ! -f $info_file ]]; then
            echo "Skipping $case_name: missing data or info beside $layer_directory" >&2
            ((skipped += 1))
            continue
        fi
        if [[ ! -f $dw_trajectory_file ]]; then
            dw_trajectory_file=""
        fi
        mkdir -p -- "$output_directory"

        job_name=$(printf '%s' "${case_name}_tamsd" | tr -c '[:alnum:]_-' '_')
        if (( ${#job_name} > 100 )); then
            job_name=${job_name:0:100}
        fi
        declare -a submit_command=(
            sbatch
            --job-name="$job_name"
            --output="$output_directory/slurm-profile-%j.out"
            --error="$output_directory/slurm-profile-%j.err"
            "$script_path"
            --worker
            "$data_file"
            "$info_file"
            "$trajectory_file"
            "$dw_trajectory_file"
            "$output_directory"
            "$time_origin_count"
            "$time_origin_stride"
            "$frame_stride"
            "$bin_width"
            "$case_name"
        )

        if (( dry_run )); then
            printf 'DRY RUN:'
            printf ' %q' "${submit_command[@]}"
            printf '\n'
        else
            "${submit_command[@]}"
        fi
        ((submitted += 1))
    done < <(find "$search_root" -type f \
        -name 'dump.layer_dynamics.*.lammpstrj' -print0 | sort -z)
done

if (( submitted == 0 )); then
    echo "No complete layer-dynamics cases were found." >&2
    exit 1
fi

echo
if (( dry_run )); then
    echo "Validated $submitted analysis job(s); nothing was submitted."
else
    echo "Submitted $submitted analysis job(s)."
fi
if (( skipped > 0 )); then
    echo "Skipped $skipped incomplete case(s)."
fi
