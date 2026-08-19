#include "pdms_filler_component.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
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

constexpr double kAtmosphereToMpa = 0.101325;
constexpr double kBondLengthAngstrom = 2.801;

struct Settings {
    std::filesystem::path data;
    std::filesystem::path source_info;
    std::filesystem::path output_directory;
    std::string mode = "auto";
    std::optional<char> direction;
    double temperature = 300.0;
    double pressure_atm = 1.0;
    double timestep_fs = 5.0;
    long long equilibration_steps = 5000000;
    double strain_rate_per_fs = 2.0e-9;
    double maximum_engineering_strain = 0.20;
    int velocity_seed = 5634;
    long long thermo_every = 1000;
    long long restart_every = 1000000;
    long long average_sample_every = 100;
    long long average_sample_count = 10;
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
    std::filesystem::path directory;
    std::string input;
    std::string submit;
    std::string info;
    std::string stress;
    std::string reference;
    std::string final_data;
};

[[noreturn]] void usage(const char* program, const std::string& message = {}) {
    if (!message.empty()) std::cerr << "Error: " << message << "\n\n";
    std::cerr
        << "Usage: " << program
        << " <data.case.npt_eq> <case.info> [options]\n\n"
        << "Options:\n"
        << "  --mode auto|parallel|all  auto: films x/y, bulk x/y/z (default)\n"
        << "  --direction x|y|z        generate only one direction\n"
        << "  --output-dir DIR         default: <source folder>/tensile\n"
        << "  --temperature K          default: 300\n"
        << "  --pressure ATM           default: 1\n"
        << "  --equilibration-steps N  default: 5000000\n"
        << "  --strain-rate RATE       engineering rate in fs^-1; default: 2e-9\n"
        << "  --max-strain VALUE       default: 0.20\n"
        << "  --velocity-seed N        default: 5634\n"
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
        throw std::runtime_error("Tensile tests require model information format version 3 or newer");
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
        throw std::runtime_error("Current PDMS tensile template requires 3 atom types");
    if (!contains_line("2 bond types"))
        throw std::runtime_error("Current PDMS tensile template requires 2 bond types");
    if (!contains_line("1 angle types") || !contains_line("1 dihedral types"))
        throw std::runtime_error("Unexpected angle or dihedral type count in data file");
}

long long parse_integer(const std::string& value, const std::string& option) {
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
        else if (option == "--mode") settings.mode = value();
        else if (option == "--direction") {
            const std::string selected = value();
            if (selected.size() != 1 || selected.find_first_not_of("xyz") != std::string::npos)
                usage(argv[0], "--direction must be x, y, or z");
            settings.direction = selected[0];
        } else if (option == "--output-dir") settings.output_directory = value();
        else if (option == "--temperature")
            settings.temperature = parse_positive(value(), option);
        else if (option == "--pressure")
            settings.pressure_atm = parse_positive(value(), option);
        else if (option == "--equilibration-steps")
            settings.equilibration_steps = parse_integer(value(), option);
        else if (option == "--strain-rate")
            settings.strain_rate_per_fs = parse_positive(value(), option);
        else if (option == "--max-strain")
            settings.maximum_engineering_strain = parse_positive(value(), option);
        else if (option == "--velocity-seed")
            settings.velocity_seed = static_cast<int>(parse_integer(value(), option));
        else if (!option.empty() && option[0] == '-')
            usage(argv[0], "unknown option " + option);
        else positional.push_back(option);
    }
    if (positional.size() != 2)
        usage(argv[0], "provide one .npt_eq file and its matching .info file");
    if (settings.mode != "auto" && settings.mode != "parallel" && settings.mode != "all")
        usage(argv[0], "--mode must be auto, parallel, or all");
    settings.data = std::filesystem::absolute(positional[0]);
    settings.source_info = std::filesystem::absolute(positional[1]);
    if (settings.output_directory.empty())
        settings.output_directory = settings.data.parent_path() / "tensile";
    else
        settings.output_directory = std::filesystem::absolute(settings.output_directory);
    return settings;
}

std::vector<char> directions(const Settings& settings, const SourceInfo& info) {
    if (settings.direction) {
        if (info.geometry == "film" && *settings.direction == 'z')
            throw std::runtime_error("Film-normal z tension is not supported by the parallel-film workflow");
        return {*settings.direction};
    }
    if (settings.mode == "parallel") return {'x', 'y'};
    if (settings.mode == "all") {
        if (info.geometry == "film")
            throw std::runtime_error("--mode all would include nonperiodic z for a film; use auto or parallel");
        return {'x', 'y', 'z'};
    }
    return info.geometry == "film" ? std::vector<char>{'x', 'y'}
                                   : std::vector<char>{'x', 'y', 'z'};
}

std::string sanitize_job_name(const std::string& value) {
    std::string result;
    for (const unsigned char c : value)
        result.push_back(std::isalnum(c) || c == '_' || c == '-' ? c : '_');
    if (result.size() > 100) result.resize(100);
    return result.empty() ? "PDMS_tensile" : result;
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

OutputFiles output_files(const Settings& settings, const SourceInfo& info, char direction) {
    OutputFiles files;
    files.directory = settings.output_directory / std::string(1, direction);
    const std::string stem = "tensile." + info.case_name + "." + direction;
    files.input = "in." + stem;
    files.submit = "submit." + stem + ".sh";
    files.info = stem + ".info";
    files.stress = "stress_strain." + info.case_name + "." + direction + ".dat";
    files.reference = "reference_box." + info.case_name + "." + direction + ".dat";
    files.final_data = "data." + info.case_name + ".tensile_" + direction + "_final";
    return files;
}

std::string relative_path(
    const std::filesystem::path& target,
    const std::filesystem::path& base
) {
    std::error_code error;
    const auto relative = std::filesystem::relative(target, base, error);
    return error ? target.string() : relative.string();
}

std::string transverse_pressure_control(char direction, const SourceInfo& info,
                                        const Settings& settings) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(6);
    if (info.geometry == "film") {
        const char transverse = direction == 'x' ? 'y' : 'x';
        output << transverse << ' ' << settings.pressure_atm << ' '
               << settings.pressure_atm << " 500.0";
        return output.str();
    }
    std::string axes;
    for (const char axis : std::string("xyz"))
        if (axis != direction) axes.push_back(axis);
    output << axes[0] << ' ' << settings.pressure_atm << ' '
           << settings.pressure_atm << " 500.0 "
           << axes[1] << ' ' << settings.pressure_atm << ' '
           << settings.pressure_atm << " 500.0 couple " << axes;
    return output.str();
}

void write_lammps_input(const Settings& settings, const SourceInfo& info,
                        char direction, const OutputFiles& files) {
    std::ofstream output(files.directory / files.input);
    if (!output)
        throw std::runtime_error("Cannot write " + (files.directory / files.input).string());

    const auto pair = pdms_filler::pair_parameters(settings.temperature);
    const double wall_cutoff = pdms_filler::repulsive_cutoff(pair);
    const long long deformation_steps = static_cast<long long>(std::ceil(
        settings.maximum_engineering_strain /
        (settings.strain_rate_per_fs * settings.timestep_fs)));
    const double realized_strain = deformation_steps * settings.strain_rate_per_fs *
                                   settings.timestep_fs;
    const std::string source_data = relative_path(settings.data, files.directory);
    const std::string direction_string(1, direction);

    output << std::fixed << std::setprecision(9)
        << "# Generated by the PDMS tensile-test generator\n"
        << "# Source network: " << info.case_name << "\n"
        << "# Loading direction: " << direction << "\n\n"
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
    pdms_filler::write_pair_matrix(output, settings.temperature, false);
    output << "\nneighbor        2 bin\n"
        << "neigh_modify    delay 5 every 1 check yes\n"
        << "timestep        " << settings.timestep_fs << "\n"
        << "thermo          " << settings.thermo_every << "\n"
        << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n\n";

    if (info.geometry == "film") {
        output << "# Fixed 300 K repulsive film walls; z is not barostatted.\n"
            << "fix             zlo_wall all wall/lj126 zlo EDGE "
            << pair.epsilon << ' ' << pair.sigma << ' ' << wall_cutoff << " units box\n"
            << "fix             zhi_wall all wall/lj126 zhi EDGE "
            << pair.epsilon << ' ' << pair.sigma << ' ' << wall_cutoff << " units box\n\n";
    }

    output << "velocity        all create " << settings.temperature << ' '
        << settings.velocity_seed << " mom yes rot yes dist gaussian\n"
        << "restart         " << settings.restart_every << " restart."
        << info.case_name << ".tensile_" << direction << ".1 restart."
        << info.case_name << ".tensile_" << direction << ".2\n\n"
        << "# Re-equilibrate the undeformed network.\n"
        << "fix             equilibrate all npt temp " << settings.temperature << ' '
        << settings.temperature << " 50.0 ";
    if (info.geometry == "film") {
        output << "x " << settings.pressure_atm << ' ' << settings.pressure_atm
               << " 500.0 y " << settings.pressure_atm << ' '
               << settings.pressure_atm << " 500.0 couple xy\n";
    } else {
        output << "iso " << settings.pressure_atm << ' ' << settings.pressure_atm
               << " 500.0\n";
    }
    output << "run             " << settings.equilibration_steps << "\n"
        << "unfix           equilibrate\n\n"
        << "# Freeze the post-equilibration box as the strain reference.\n"
        << "run             0\n"
        << "variable        Lx0 equal $(lx:%.12g)\n"
        << "variable        Ly0 equal $(ly:%.12g)\n"
        << "variable        Lz0 equal $(lz:%.12g)\n"
        << "print           \"Lx0_A Ly0_A Lz0_A\" file " << files.reference
        << " screen no\n"
        << "print           \"${Lx0} ${Ly0} ${Lz0}\" append " << files.reference
        << " screen no\n"
        << "reset_timestep  0\n\n"
        << "variable        time_ns equal step*dt*1.0e-6\n"
        << "variable        lambda_x equal lx/v_Lx0\n"
        << "variable        lambda_y equal ly/v_Ly0\n"
        << "variable        lambda_z equal lz/v_Lz0\n"
        << "variable        strain_x equal v_lambda_x-1.0\n"
        << "variable        strain_y equal v_lambda_y-1.0\n"
        << "variable        strain_z equal v_lambda_z-1.0\n"
        << "variable        sample_temperature equal temp\n"
        << "variable        pressure_x equal pxx\n"
        << "variable        pressure_y equal pyy\n"
        << "variable        pressure_z equal pzz\n";

    if (info.geometry == "film") {
        output << "variable        material_thickness equal " << info.film_thickness << "\n"
            << "variable        material_volume equal lx*ly*v_material_thickness\n"
            << "variable        volume_correction equal vol/v_material_volume\n";
    } else {
        output << "variable        material_volume equal vol\n"
            << "variable        volume_correction equal 1.0\n";
    }
    output << "variable        sigma_x equal -pxx*" << kAtmosphereToMpa
        << "*v_volume_correction\n"
        << "variable        sigma_y equal -pyy*" << kAtmosphereToMpa
        << "*v_volume_correction\n"
        << "variable        sigma_z equal -pzz*" << kAtmosphereToMpa
        << "*v_volume_correction\n";

    if (info.geometry == "film" && direction == 'x') {
        output << "variable        reduced_strain equal v_lambda_x^2-v_lambda_y^2\n"
            << "variable        stress_difference equal v_sigma_x-v_sigma_y\n";
    } else if (info.geometry == "film") {
        output << "variable        reduced_strain equal v_lambda_y^2-v_lambda_x^2\n"
            << "variable        stress_difference equal v_sigma_y-v_sigma_x\n";
    } else {
        std::string transverse;
        for (const char axis : std::string("xyz"))
            if (axis != direction) transverse.push_back(axis);
        output << "variable        reduced_strain equal v_lambda_" << direction
            << "^2-0.5*(v_lambda_" << transverse[0] << "^2+v_lambda_"
            << transverse[1] << "^2)\n"
            << "variable        stress_difference equal v_sigma_" << direction
            << "-0.5*(v_sigma_" << transverse[0] << "+v_sigma_"
            << transverse[1] << ")\n";
    }

    output << "\nthermo_style    custom step v_time_ns temp lx ly lz"
        << " v_lambda_x v_lambda_y v_lambda_z pxx pyy pzz"
        << " v_sigma_x v_sigma_y v_sigma_z v_reduced_strain v_stress_difference\n\n"
        << "# The table contains block-averaged values every "
        << settings.average_sample_every * settings.average_sample_count
        << " steps.\n"
        << "fix             samples all ave/time " << settings.average_sample_every
        << ' ' << settings.average_sample_count << ' '
        << settings.average_sample_every * settings.average_sample_count
        << " v_time_ns v_lambda_x v_lambda_y v_lambda_z"
        << " v_strain_x v_strain_y v_strain_z v_sample_temperature"
        << " v_pressure_x v_pressure_y v_pressure_z"
        << " v_sigma_x v_sigma_y v_sigma_z v_material_volume"
        << " v_volume_correction v_reduced_strain v_stress_difference"
        << " file " << files.stress << " ave one format \" %.16g\"\n\n"
        << "# Uniaxial deformation; only transverse periodic dimensions are barostatted.\n"
        << "fix             integrate all npt temp " << settings.temperature << ' '
        << settings.temperature << " 50.0 "
        << transverse_pressure_control(direction, info, settings) << "\n"
        << "fix             deform_box all deform 1 " << direction << " erate "
        << std::scientific << std::setprecision(12) << settings.strain_rate_per_fs
        << std::fixed << std::setprecision(9) << " units box remap x\n"
        << "run             " << deformation_steps << "\n\n"
        << "unfix           deform_box\n"
        << "unfix           integrate\n"
        << "unfix           samples\n"
        << "write_data      " << files.final_data << " nocoeff\n"
        << "print           \"Tensile test completed: direction " << direction
        << ", target engineering strain " << settings.maximum_engineering_strain
        << ", realized " << realized_strain << "\"\n";
}

void write_submit(const SourceInfo& info, char direction, const OutputFiles& files) {
    std::ofstream output(files.directory / files.submit);
    if (!output)
        throw std::runtime_error("Cannot write " + (files.directory / files.submit).string());
    output << "#!/bin/bash\n"
        << "#SBATCH --job-name="
        << sanitize_job_name(info.case_name + "_ten_" + direction) << "\n"
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
        << "OUTPUT=" << shell_single_quote("out.tensile." + info.case_name + "." + direction)
        << "\n"
        << "srun lmp -in \"$INPUT\" > \"$OUTPUT\"\n";
}

void write_info_file(const Settings& settings, const SourceInfo& info,
                     char direction, const OutputFiles& files) {
    std::ofstream output(files.directory / files.info);
    if (!output)
        throw std::runtime_error("Cannot write " + (files.directory / files.info).string());
    const long long deformation_steps = static_cast<long long>(std::ceil(
        settings.maximum_engineering_strain /
        (settings.strain_rate_per_fs * settings.timestep_fs)));
    const double realized_strain = deformation_steps * settings.strain_rate_per_fs *
                                   settings.timestep_fs;
    output << std::fixed << std::setprecision(10)
        << "{\n"
        << "  \"format\": \"pdms-elastomer-tensile-info\",\n"
        << "  \"format_version\": 1,\n"
        << "  \"case_name\": \"" << json_escape(info.case_name) << "\",\n"
        << "  \"architecture\": \"" << json_escape(info.architecture) << "\",\n"
        << "  \"geometry\": \"" << info.geometry << "\",\n"
        << "  \"loading_direction\": \"" << direction << "\",\n"
        << "  \"source\": {\n"
        << "    \"data\": \"" << json_escape(relative_path(settings.data, files.directory))
        << "\",\n"
        << "    \"model_info\": \""
        << json_escape(relative_path(settings.source_info, files.directory)) << "\",\n"
        << "    \"model_info_format_version\": " << info.format_version << "\n"
        << "  },\n"
        << "  \"files\": {\n"
        << "    \"lammps_input\": \"" << files.input << "\",\n"
        << "    \"slurm_submit\": \"" << files.submit << "\",\n"
        << "    \"stress_strain\": \"" << files.stress << "\",\n"
        << "    \"runtime_reference_box\": \"" << files.reference << "\",\n"
        << "    \"final_data\": \"" << files.final_data << "\"\n"
        << "  },\n"
        << "  \"protocol\": {\n"
        << "    \"temperature_K\": " << settings.temperature << ",\n"
        << "    \"transverse_pressure_atm\": " << settings.pressure_atm << ",\n"
        << "    \"timestep_fs\": " << settings.timestep_fs << ",\n"
        << "    \"equilibration_steps\": " << settings.equilibration_steps << ",\n"
        << "    \"equilibration_duration_ns\": "
        << settings.equilibration_steps * settings.timestep_fs * 1.0e-6 << ",\n"
        << "    \"engineering_strain_rate_per_fs\": "
        << settings.strain_rate_per_fs << ",\n"
        << "    \"requested_maximum_engineering_strain\": "
        << settings.maximum_engineering_strain << ",\n"
        << "    \"deformation_steps\": " << deformation_steps << ",\n"
        << "    \"deformation_duration_ns\": "
        << deformation_steps * settings.timestep_fs * 1.0e-6 << ",\n"
        << "    \"realized_maximum_engineering_strain\": " << realized_strain << ",\n"
        << "    \"velocity_seed\": " << settings.velocity_seed << ",\n"
        << "    \"trajectory_dump_enabled\": false,\n"
        << "    \"reference_box_definition\": "
        << "\"captured after tensile equilibration and before reset_timestep\"\n"
        << "  },\n"
        << "  \"film\": {\n"
        << "    \"nominal_material_thickness_angstrom\": ";
    if (info.geometry == "film") output << info.film_thickness;
    else output << "null";
    output << ",\n    \"wall_cutoff_per_side_angstrom\": ";
    if (info.geometry == "film") output << info.wall_cutoff;
    else output << "null";
    output << ",\n    \"wall_separation_fixed_during_deformation\": "
        << (info.geometry == "film" ? "true" : "null") << ",\n"
        << "    \"material_volume_definition\": \""
        << (info.geometry == "film" ? "Lx*Ly*nominal_material_thickness" : "box volume")
        << "\",\n"
        << "    \"stress_volume_correction\": \""
        << (info.geometry == "film" ? "box volume / material volume = Lz/H" : "1")
        << "\"\n"
        << "  },\n"
        << "  \"analysis_contract\": {\n"
        << "    \"pressure_sign\": \"LAMMPS pressure is positive in compression\",\n"
        << "    \"pressure_conversion_MPa_per_atm\": " << kAtmosphereToMpa << ",\n"
        << "    \"primary_fit\": \"stress_difference = intercept + G*reduced_strain\",\n"
        << "    \"default_fit_engineering_strain_range\": [0.0, 0.05],\n"
        << "    \"stress_table_columns\": [\"time_ns\", \"lambda_x\", \"lambda_y\","
        << " \"lambda_z\", \"strain_x\", \"strain_y\", \"strain_z\","
        << " \"temperature_K\", \"Pxx_atm\", \"Pyy_atm\", \"Pzz_atm\","
        << " \"sigma_x_MPa\", \"sigma_y_MPa\", \"sigma_z_MPa\","
        << " \"material_volume_A3\", \"volume_correction\","
        << " \"reduced_strain\", \"stress_difference_MPa\"]\n"
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

        for (const char direction : directions(settings, info)) {
            const OutputFiles files = output_files(settings, info, direction);
            std::filesystem::create_directories(files.directory, error);
            if (error)
                throw std::runtime_error(
                    "Cannot create output directory " + files.directory.string() +
                    ": " + error.message());
            write_lammps_input(settings, info, direction, files);
            write_submit(info, direction, files);
            write_info_file(settings, info, direction, files);
            std::filesystem::permissions(
                files.directory / files.submit,
                std::filesystem::perms::owner_exec |
                    std::filesystem::perms::group_exec |
                    std::filesystem::perms::others_exec,
                std::filesystem::perm_options::add,
                error);
            std::cout << "Generated " << direction << " tensile test in "
                      << files.directory << "\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "tensile_test_generator: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
