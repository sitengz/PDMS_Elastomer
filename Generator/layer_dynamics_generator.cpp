#include "pdms_filler_component.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kTemperatureK = 300.0;
constexpr double kTimestepFs = 5.0;
constexpr double kBondLengthAngstrom = 2.801;

struct Settings {
    std::filesystem::path data;
    std::filesystem::path source_info;
    std::filesystem::path output_directory;
    long long surface_equilibration_steps = 1000000;
    long long debye_waller_steps = 20000;
    long long debye_waller_dump_every = 20;
    long long production_steps = 5000000;
    long long early_steps = 1000000;
    long long early_dump_every = 1000;
    long long long_dump_every = 5000;
    long long thermo_every = 1000;
    long long restart_every = 1000000;
    double film_vacuum_padding_per_side = 50.0;
};

struct SourceInfo {
    std::string case_name;
    std::string geometry;
    std::string architecture;
    int format_version = 0;
    double film_thickness = 0.0;
    double wall_cutoff = 0.0;
};

struct OutputFiles {
    std::string input;
    std::string submit;
    std::string info;
    std::string debye_waller_trajectory;
    std::string msd_trajectory;
    std::string equilibrated_data;
    std::string final_data;
};

[[noreturn]] void usage(const char* program, const std::string& message = {}) {
    if (!message.empty()) std::cerr << "Error: " << message << "\n\n";
    std::cerr
        << "Usage: " << program
        << " <data.case.npt_eq> <case.info> [options]\n\n"
        << "Options:\n"
        << "  --output-dir DIR          default: <source folder>/layer_dynamics\n"
        << "  --surface-equilibration-steps N  default: 1000000 (5 ns)\n"
        << "  --dw-steps N              default: 20000 (100 ps)\n"
        << "  --dw-dump-every N         default: 20 steps (0.1 ps)\n"
        << "  --production-steps N      default: 5000000 (25 ns)\n"
        << "  --early-steps N           default: 1000000 (5 ns)\n"
        << "  --early-dump-every N      default: 1000 steps (5 ps)\n"
        << "  --long-dump-every N       default: 5000 steps (25 ps)\n"
        << "  --film-padding X          vacuum added to each z side; default: 50 A\n"
        << "  --thermo-every N          default: 1000\n"
        << "  --restart-every N         default: 1000000\n"
        << "  --help\n";
    std::exit(message.empty() ? 0 : 2);
}

std::string read_text(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read file: " + path.string());
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::string regex_escape(const std::string& text) {
    static const std::regex special(R"([.^$|()\[\]{}*+?\\])");
    return std::regex_replace(text, special, R"(\$&)");
}

std::optional<std::string> json_string(
    const std::string& text,
    const std::string& key
) {
    const std::regex pattern(
        "\\\"" + regex_escape(key) + "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) return std::nullopt;
    return match[1].str();
}

std::optional<double> json_number(
    const std::string& text,
    const std::string& key
) {
    const std::regex pattern(
        "\\\"" + regex_escape(key) +
        "\\\"\\s*:\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) return std::nullopt;
    return std::stod(match[1].str());
}

SourceInfo parse_source_info(const std::filesystem::path& path) {
    const std::string text = read_text(path);
    const auto format = json_string(text, "format");
    if (!format || *format != "pdms-elastomer-model-info")
        throw std::runtime_error("Unsupported model information format in " + path.string());

    SourceInfo info;
    info.case_name = json_string(text, "case_name").value_or("");
    info.geometry = json_string(text, "geometry").value_or("");
    info.architecture = json_string(text, "strand_topology").value_or("unknown");
    info.format_version = static_cast<int>(json_number(text, "format_version").value_or(0));
    info.film_thickness = json_number(text, "film_thickness_angstrom").value_or(0.0);
    info.wall_cutoff =
        json_number(text, "film_wall_cutoff_per_side_angstrom").value_or(0.0);

    if (info.format_version < 3)
        throw std::runtime_error(
            "Layer dynamics requires model information format version 3 or newer");
    if (info.case_name.empty())
        throw std::runtime_error("Missing case_name in " + path.string());
    if (info.geometry != "bulk" && info.geometry != "film")
        throw std::runtime_error("geometry must be bulk or film in " + path.string());
    if (info.geometry == "film" &&
        (!(info.film_thickness > 0.0) || !(info.wall_cutoff > 0.0)))
        throw std::runtime_error("Film thickness and wall cutoff must be positive");
    return info;
}

void validate_data(const std::filesystem::path& path, const SourceInfo& info) {
    if (!std::filesystem::is_regular_file(path))
        throw std::runtime_error("Missing equilibrated data file: " + path.string());
    const std::string expected = "data." + info.case_name + ".npt_eq";
    if (path.filename() != expected)
        throw std::runtime_error(
            "Expected equilibrated data filename " + expected + ", received " +
            path.filename().string());

    const std::string text = read_text(path);
    const auto contains_line = [&text](const std::string& value) {
        return text.find("\n" + value + "\n") != std::string::npos ||
               text.rfind(value + "\n", 0) == 0;
    };
    if (!contains_line("3 atom types"))
        throw std::runtime_error("Current PDMS dynamics template requires 3 atom types");
    if (!contains_line("2 bond types"))
        throw std::runtime_error("Current PDMS dynamics template requires 2 bond types");
    if (!contains_line("1 angle types") || !contains_line("1 dihedral types"))
        throw std::runtime_error("Unexpected angle or dihedral type count in data file");
    if (text.find("\nVelocities\n") == std::string::npos)
        throw std::runtime_error(
            "The equilibrated data file has no Velocities section; layer dynamics "
            "preserves the 300 K velocities from the source state");
}

long long parse_positive_integer(const std::string& value, const std::string& option) {
    std::size_t used = 0;
    const long long result = std::stoll(value, &used);
    if (used != value.size() || result <= 0)
        throw std::runtime_error(option + " must be a positive integer");
    return result;
}

double parse_positive(const std::string& value, const std::string& option) {
    std::size_t used = 0;
    const double result = std::stod(value, &used);
    if (used != value.size() || !(result > 0.0) || !std::isfinite(result))
        throw std::runtime_error(option + " must be a positive finite number");
    return result;
}

Settings parse_arguments(int argc, char** argv) {
    if (argc == 1) usage(argv[0], "missing input files");
    Settings settings;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        const auto value = [&]() -> std::string {
            if (++i >= argc) usage(argv[0], "missing value for " + option);
            return argv[i];
        };
        if (option == "--help" || option == "-h") usage(argv[0]);
        else if (option == "--output-dir") settings.output_directory = value();
        else if (option == "--surface-equilibration-steps")
            settings.surface_equilibration_steps = parse_positive_integer(value(), option);
        else if (option == "--dw-steps")
            settings.debye_waller_steps = parse_positive_integer(value(), option);
        else if (option == "--dw-dump-every")
            settings.debye_waller_dump_every = parse_positive_integer(value(), option);
        else if (option == "--production-steps")
            settings.production_steps = parse_positive_integer(value(), option);
        else if (option == "--early-steps")
            settings.early_steps = parse_positive_integer(value(), option);
        else if (option == "--early-dump-every")
            settings.early_dump_every = parse_positive_integer(value(), option);
        else if (option == "--long-dump-every")
            settings.long_dump_every = parse_positive_integer(value(), option);
        else if (option == "--film-padding")
            settings.film_vacuum_padding_per_side = parse_positive(value(), option);
        else if (option == "--thermo-every")
            settings.thermo_every = parse_positive_integer(value(), option);
        else if (option == "--restart-every")
            settings.restart_every = parse_positive_integer(value(), option);
        else if (!option.empty() && option[0] == '-')
            usage(argv[0], "unknown option " + option);
        else positional.push_back(option);
    }
    if (positional.size() != 2)
        usage(argv[0], "provide one .npt_eq file and its matching .info file");
    if (settings.early_steps > settings.production_steps)
        usage(argv[0], "--early-steps cannot exceed --production-steps");
    if (settings.debye_waller_steps % settings.debye_waller_dump_every != 0)
        usage(argv[0], "--dw-steps must be divisible by --dw-dump-every");
    if (settings.early_steps % settings.early_dump_every != 0)
        usage(argv[0], "--early-steps must be divisible by --early-dump-every");
    if (settings.early_steps % settings.long_dump_every != 0 ||
        settings.production_steps % settings.long_dump_every != 0)
        usage(argv[0], "early and production steps must be divisible by --long-dump-every");

    settings.data = std::filesystem::absolute(positional[0]);
    settings.source_info = std::filesystem::absolute(positional[1]);
    if (settings.output_directory.empty())
        settings.output_directory = settings.data.parent_path() / "layer_dynamics";
    else
        settings.output_directory = std::filesystem::absolute(settings.output_directory);
    return settings;
}

std::string sanitize_job_name(const std::string& value) {
    std::string result;
    for (const unsigned char c : value)
        result.push_back(std::isalnum(c) || c == '_' || c == '-' ? c : '_');
    if (result.size() > 100) result.resize(100);
    return result.empty() ? "PDMS_layer_dyn" : result;
}

std::string shell_single_quote(const std::string& value) {
    std::string result = "'";
    for (const char c : value) result += c == '\'' ? "'\\''" : std::string(1, c);
    return result + "'";
}

std::string json_escape(const std::string& value) {
    std::ostringstream output;
    for (const unsigned char c : value) {
        if (c == '"') output << "\\\"";
        else if (c == '\\') output << "\\\\";
        else if (c == '\n') output << "\\n";
        else if (c == '\r') output << "\\r";
        else if (c == '\t') output << "\\t";
        else if (c < 0x20)
            output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                   << static_cast<int>(c) << std::dec << std::setfill(' ');
        else output << static_cast<char>(c);
    }
    return output.str();
}

std::string relative_path(
    const std::filesystem::path& target,
    const std::filesystem::path& base
) {
    std::error_code error;
    const auto relative = std::filesystem::relative(target, base, error);
    return error ? target.string() : relative.string();
}

OutputFiles output_files(const SourceInfo& info) {
    OutputFiles files;
    files.input = "in.layer_dynamics." + info.case_name;
    files.submit = "submit.layer_dynamics." + info.case_name + ".sh";
    files.info = "layer_dynamics." + info.case_name + ".info";
    files.debye_waller_trajectory =
        "dump.debye_waller." + info.case_name + ".lammpstrj";
    files.msd_trajectory =
        "dump.layer_dynamics." + info.case_name + ".lammpstrj";
    files.equilibrated_data = info.geometry == "film"
        ? "data." + info.case_name + ".free_surface_eq"
        : "data." + info.case_name + ".layer_dynamics_eq";
    files.final_data = "data." + info.case_name + ".layer_dynamics_final";
    return files;
}

void write_lammps_input(const Settings& settings, const SourceInfo& info,
                        const OutputFiles& files) {
    std::ofstream output(settings.output_directory / files.input);
    if (!output)
        throw std::runtime_error(
            "Cannot write " + (settings.output_directory / files.input).string());

    const long long long_steps = settings.production_steps - settings.early_steps;
    const std::string source_data = relative_path(settings.data, settings.output_directory);
    const pdms_filler::PairParameters guard_wall =
        pdms_filler::pair_parameters(kTemperatureK);
    const double guard_wall_cutoff = pdms_filler::repulsive_cutoff(guard_wall);

    output << std::fixed << std::setprecision(9)
        << "# Generated by the PDMS layer-dynamics generator\n"
        << "# Source network: " << info.case_name << "\n"
        << "# 300 K NVT trajectory with inherited equilibrated velocities\n\n"
        << "units           real\n"
        << "boundary        p p " << (info.geometry == "film" ? "f" : "p") << "\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     harmonic\n"
        << "dihedral_style  nharmonic\n"
        << "special_bonds   lj 0 0 0.5\n"
        << "pair_style      lj/gromacs 12 15\n"
        << "comm_modify     cutoff 25\n"
        << "read_data       " << source_data << "\n\n"
        << "mass            1 74.000000000\n"
        << "mass            2 74.000000000\n"
        << "mass            3 74.000000000\n\n"
        << "bond_coeff      1 115.4086 " << kBondLengthAngstrom << "\n"
        << "bond_coeff      2 115.4086 " << kBondLengthAngstrom << "\n"
        << "angle_coeff     1 64.62431 111.623\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n\n";
    pdms_filler::write_pair_matrix(output, kTemperatureK, false);
    output << "\nneighbor        2 bin\n"
        << "neigh_modify    delay 5 every 1 check yes\n"
        << "timestep        " << kTimestepFs << "\n";

    output << "\n# The .npt_eq Velocities section is retained; no velocity reinitialization.\n";
    if (info.geometry == "film") {
        output << "# Add vacuum on both z sides, then place remote repulsive guard walls\n"
            << "# at the expanded edges. The original material-adjacent walls are not restored.\n"
            << "change_box      all z delta -" << settings.film_vacuum_padding_per_side
            << ' ' << settings.film_vacuum_padding_per_side << " units box\n"
            << "fix             zlo_guard all wall/lj126 zlo EDGE "
            << guard_wall.epsilon << ' ' << guard_wall.sigma << ' '
            << guard_wall_cutoff << " units box\n"
            << "fix             zhi_guard all wall/lj126 zhi EDGE "
            << guard_wall.epsilon << ' ' << guard_wall.sigma << ' '
            << guard_wall_cutoff << " units box\n"
            << "fix_modify      zlo_guard energy yes\n"
            << "fix_modify      zhi_guard energy yes\n";
    }
    output << "\nthermo          " << settings.thermo_every << "\n"
        << "thermo_style    custom step time temp density lx ly lz"
        << " etotal epair ebond eangle edihed\n"
        << "thermo_modify   flush yes\n\n"
        << "restart         " << settings.restart_every << " restart."
        << info.case_name << ".layer_dynamics.1 restart."
        << info.case_name << ".layer_dynamics.2\n\n"
        << "# Equilibrate after expansion and remote guard-wall placement.\n"
        << "fix             equilibrate all nvt temp " << kTemperatureK << ' '
        << kTemperatureK << " 50.0\n"
        << "run             " << settings.surface_equilibration_steps << "\n"
        << "unfix           equilibrate\n"
        << "write_data      " << files.equilibrated_data << " nocoeff\n\n"
        << "# Short, high-frequency trajectory for the Debye-Waller displacement.\n"
        << "reset_timestep  0\n"
        << "compute         global_dw all msd com yes\n"
        << "thermo          " << settings.thermo_every << "\n"
        << "thermo_style    custom step time temp density lx ly lz"
        << " c_global_dw[1] c_global_dw[2] c_global_dw[3] c_global_dw[4]"
        << " etotal epair ebond eangle edihed\n"
        << "thermo_modify   flush yes\n\n"
        << "# First frame is the Debye-Waller displacement origin.\n"
        << "dump            dw all custom " << settings.debye_waller_dump_every << ' '
        << files.debye_waller_trajectory << " id mol type x y z ix iy iz\n"
        << "dump_modify     dw first yes sort id\n"
        << "fix             dw_integrate all nvt temp " << kTemperatureK << ' '
        << kTemperatureK << " 50.0\n\n"
        << "run             " << settings.debye_waller_steps << "\n"
        << "unfix           dw_integrate\n"
        << "undump          dw\n"
        << "uncompute       global_dw\n\n"
        << "# Reset the reference for the independent long-time layer-MSD trajectory.\n"
        << "reset_timestep  0\n"
        << "compute         global_msd all msd com yes\n"
        << "thermo_style    custom step time temp density lx ly lz"
        << " c_global_msd[1] c_global_msd[2] c_global_msd[3] c_global_msd[4]"
        << " etotal epair ebond eangle edihed\n"
        << "dump            msd all custom " << settings.early_dump_every << ' '
        << files.msd_trajectory << " id mol type x y z ix iy iz\n"
        << "dump_modify     msd first yes sort id\n"
        << "fix             msd_integrate all nvt temp " << kTemperatureK << ' '
        << kTemperatureK << " 50.0\n\n"
        << "# MSD early-time window: dense sampling through "
        << settings.early_steps * kTimestepFs * 1.0e-6 << " ns.\n"
        << "run             " << settings.early_steps << "\n";
    if (long_steps > 0) {
        output << "\n# Production long-time window: coarser sampling through "
            << settings.production_steps * kTimestepFs * 1.0e-6 << " ns.\n"
            << "dump_modify     msd every " << settings.long_dump_every
            << " first no\n"
            << "run             " << long_steps << "\n";
    }
    output << "\nunfix           msd_integrate\n"
        << "undump          msd\n"
        << "# Remove compute-dependent thermo fields before deleting the MSD compute.\n"
        << "thermo_style    custom step time temp density lx ly lz"
        << " etotal epair ebond eangle edihed\n"
        << "thermo_modify   flush yes\n"
        << "uncompute       global_msd\n"
        << "write_data      " << files.final_data << " nocoeff\n"
        << "print           \"Layer-dynamics production completed: "
        << settings.production_steps << " steps, "
        << settings.production_steps * kTimestepFs * 1.0e-6 << " ns\"\n";
}

void write_submit(const SourceInfo& info, const OutputFiles& files,
                  const std::filesystem::path& directory) {
    std::ofstream output(directory / files.submit);
    if (!output)
        throw std::runtime_error("Cannot write " + (directory / files.submit).string());
    output << "#!/bin/bash\n"
        << "#SBATCH --job-name=" << sanitize_job_name(info.case_name + "_layer_dyn") << "\n"
        << "#SBATCH --time=48:00:00\n"
        << "#SBATCH --nodes=1\n"
        << "#SBATCH --ntasks-per-node=96\n"
        << "#SBATCH --mem=200G\n"
        << "#SBATCH --partition=nova\n"
        << "#SBATCH --mail-user=siteng@iastate.edu\n"
        << "#SBATCH --mail-type=END,FAIL\n"
        << "#SBATCH --output=slurm-%j.out\n"
        << "#SBATCH --error=slurm-%j.err\n\n"
        << "set -euo pipefail\n"
        << "cd -- \"${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}\"\n\n"
        << "module purge\n"
        << "module load intel/22.3.1\n"
        << "module load mpi/2021.7.1\n"
        << "module load lammps/20230802.2-py310-openmpi4-ezoqd7f\n\n"
        << "export OMP_NUM_THREADS=1\n\n"
        << "INPUT=" << shell_single_quote(files.input) << "\n"
        << "OUTPUT=" << shell_single_quote("out.layer_dynamics." + info.case_name) << "\n"
        << "srun lmp -in \"$INPUT\" > \"$OUTPUT\"\n";
}

void write_info_file(const Settings& settings, const SourceInfo& info,
                     const OutputFiles& files) {
    std::ofstream output(settings.output_directory / files.info);
    if (!output)
        throw std::runtime_error(
            "Cannot write " + (settings.output_directory / files.info).string());
    const long long expected_frames =
        settings.early_steps / settings.early_dump_every + 1 +
        (settings.production_steps - settings.early_steps) / settings.long_dump_every;
    const long long expected_dw_frames =
        settings.debye_waller_steps / settings.debye_waller_dump_every + 1;
    const pdms_filler::PairParameters guard_wall =
        pdms_filler::pair_parameters(kTemperatureK);
    const double guard_wall_cutoff = pdms_filler::repulsive_cutoff(guard_wall);
    output << std::fixed << std::setprecision(10)
        << "{\n"
        << "  \"format\": \"pdms-elastomer-layer-dynamics-info\",\n"
        << "  \"format_version\": 2,\n"
        << "  \"case_name\": \"" << json_escape(info.case_name) << "\",\n"
        << "  \"architecture\": \"" << json_escape(info.architecture) << "\",\n"
        << "  \"geometry\": \"" << info.geometry << "\",\n"
        << "  \"source\": {\n"
        << "    \"data\": \""
        << json_escape(relative_path(settings.data, settings.output_directory)) << "\",\n"
        << "    \"model_info\": \""
        << json_escape(relative_path(settings.source_info, settings.output_directory)) << "\",\n"
        << "    \"model_info_format_version\": " << info.format_version << "\n"
        << "  },\n"
        << "  \"files\": {\n"
        << "    \"lammps_input\": \"" << files.input << "\",\n"
        << "    \"slurm_submit\": \"" << files.submit << "\",\n"
        << "    \"debye_waller_trajectory\": \""
        << files.debye_waller_trajectory << "\",\n"
        << "    \"msd_trajectory\": \"" << files.msd_trajectory << "\",\n"
        << "    \"equilibrated_data\": \"" << files.equilibrated_data << "\",\n"
        << "    \"final_data\": \"" << files.final_data << "\"\n"
        << "  },\n"
        << "  \"protocol\": {\n"
        << "    \"ensemble\": \"NVT\",\n"
        << "    \"temperature_K\": " << kTemperatureK << ",\n"
        << "    \"temperature_damping_fs\": 50.0,\n"
        << "    \"timestep_fs\": " << kTimestepFs << ",\n"
        << "    \"surface_equilibration_steps\": "
        << settings.surface_equilibration_steps << ",\n"
        << "    \"surface_equilibration_duration_ns\": "
        << settings.surface_equilibration_steps * kTimestepFs * 1.0e-6 << ",\n"
        << "    \"debye_waller_steps\": " << settings.debye_waller_steps << ",\n"
        << "    \"debye_waller_duration_ps\": "
        << settings.debye_waller_steps * kTimestepFs * 1.0e-3 << ",\n"
        << "    \"debye_waller_dump_every_steps\": "
        << settings.debye_waller_dump_every << ",\n"
        << "    \"debye_waller_dump_every_ps\": "
        << settings.debye_waller_dump_every * kTimestepFs * 1.0e-3 << ",\n"
        << "    \"expected_debye_waller_frames\": " << expected_dw_frames << ",\n"
        << "    \"production_steps\": " << settings.production_steps << ",\n"
        << "    \"production_duration_ns\": "
        << settings.production_steps * kTimestepFs * 1.0e-6 << ",\n"
        << "    \"source_velocities_retained\": true,\n"
        << "    \"velocity_reinitialization\": false,\n"
        << "    \"box_dimensions_fixed_during_production\": true,\n"
        << "    \"trajectory_coordinates\": \"wrapped x y z plus ix iy iz\",\n"
        << "    \"each_trajectory_first_frame_is_origin\": true,\n"
        << "    \"msd_trajectory_schedule\": [\n"
        << "      {\"start_step\": 0, \"end_step\": " << settings.early_steps
        << ", \"dump_every_steps\": " << settings.early_dump_every
        << ", \"dump_every_ps\": "
        << settings.early_dump_every * kTimestepFs * 1.0e-3 << "},\n"
        << "      {\"start_step\": " << settings.early_steps
        << ", \"end_step\": " << settings.production_steps
        << ", \"dump_every_steps\": " << settings.long_dump_every
        << ", \"dump_every_ps\": "
        << settings.long_dump_every * kTimestepFs * 1.0e-3 << "}\n"
        << "    ],\n"
        << "    \"expected_msd_trajectory_frames\": " << expected_frames << "\n"
        << "  },\n"
        << "  \"film\": {\n"
        << "    \"nominal_material_thickness_angstrom\": ";
    if (info.geometry == "film") output << info.film_thickness;
    else output << "null";
    output << ",\n    \"wall_cutoff_per_side_angstrom\": ";
    if (info.geometry == "film") output << info.wall_cutoff;
    else output << "null";
    output << ",\n    \"source_walls_recreated\": "
        << (info.geometry == "film" ? "false" : "null") << ",\n"
        << "    \"remote_guard_walls\": "
        << (info.geometry == "film" ? "true" : "null") << ",\n"
        << "    \"guard_walls_at_expanded_box_edges\": "
        << (info.geometry == "film" ? "true" : "null") << ",\n"
        << "    \"guard_wall_style\": ";
    if (info.geometry == "film") output << "\"wall/lj126\"";
    else output << "null";
    output << ",\n    \"guard_wall_epsilon_kcal_per_mol\": ";
    if (info.geometry == "film") output << guard_wall.epsilon;
    else output << "null";
    output << ",\n    \"guard_wall_sigma_angstrom\": ";
    if (info.geometry == "film") output << guard_wall.sigma;
    else output << "null";
    output << ",\n    \"guard_wall_cutoff_angstrom\": ";
    if (info.geometry == "film") output << guard_wall_cutoff;
    else output << "null";
    output << ",\n"
        << "    \"free_surfaces\": "
        << (info.geometry == "film" ? "true" : "null") << ",\n"
        << "    \"vacuum_padding_per_side_angstrom\": ";
    if (info.geometry == "film") output << settings.film_vacuum_padding_per_side;
    else output << "null";
    output << "\n"
        << "  },\n"
        << "  \"analysis_contract\": {\n"
        << "    \"analyzer\": \"network_profile_analyzer\",\n"
        << "    \"debye_waller_trajectory_option\": \"--dw-trajectory "
        << files.debye_waller_trajectory << "\",\n"
        << "    \"msd_trajectory_option\": \"--trajectory "
        << files.msd_trajectory << "\",\n"
        << "    \"layer_assignment\": \"component-1 beads grouped by first-frame z\",\n"
        << "    \"displacement_reference\": \"first trajectory frame\",\n"
        << "    \"drift_correction\": \"whole-system center-of-mass displacement\",\n"
        << "    \"reported_components\": [\"u2_xy\", \"u2_3D\","
        << " \"local_stiffness_xy\", \"local_stiffness_3D\","
        << " \"MSD_x\", \"MSD_y\", \"MSD_z\", \"MSD_parallel\","
        << " \"MSD_total\", \"D_xy\", \"D_3D\"],\n"
        << "    \"diffusion_definitions\": {\"D_xy\": \"slope(MSD_parallel)/4\","
        << " \"D_3D\": \"slope(MSD_total)/6\"}\n"
        << "  }\n"
        << "}\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Settings settings = parse_arguments(argc, argv);
        const SourceInfo info = parse_source_info(settings.source_info);
        validate_data(settings.data, info);

        std::error_code error;
        std::filesystem::create_directories(settings.output_directory, error);
        if (error)
            throw std::runtime_error(
                "Cannot create output directory " + settings.output_directory.string() +
                ": " + error.message());

        const OutputFiles files = output_files(info);
        write_lammps_input(settings, info, files);
        write_submit(info, files, settings.output_directory);
        write_info_file(settings, info, files);
        std::filesystem::permissions(
            settings.output_directory / files.submit,
            std::filesystem::perms::owner_exec |
                std::filesystem::perms::group_exec |
                std::filesystem::perms::others_exec,
            std::filesystem::perm_options::add,
            error);

        std::cout << "Generated layer-dynamics test in "
                  << settings.output_directory << "\n";
        std::cout << "  " << files.input << "\n"
                  << "  " << files.info << "\n"
                  << "  " << files.submit << "\n";
    } catch (const std::exception& error) {
        std::cerr << "layer_dynamics_generator: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
