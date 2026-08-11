#include "network_common.hpp"

#include <iostream>

namespace {
using namespace pdms_analysis;

struct Options {
    std::string data_file;
    std::string info_file;
    std::string trajectory_file;
    std::string output_directory;
    double bin_width = 5.0;
    long long frame_stride = 1;
};

struct BinCounts {
    long long crosslink_bonds = 0;
    long long strand_crosslink_bonds = 0;
    long long moderator_crosslink_bonds = 0;
    long long junctions = 0;
    long long active_strands = 0;
    long long dangling_ends = 0;
    long long dangling_loops = 0;
    long long self_loops = 0;
    long long isolated_parents = 0;
};

struct DumpFrame {
    long long timestep = 0;
    Box box;
    std::vector<Vec3> unwrapped;
};

std::vector<std::string> words(const std::string &line) {
    std::vector<std::string> result;
    std::istringstream fields(line);
    std::string value;
    while (fields >> value) result.push_back(value);
    return result;
}

class DumpReader {
  public:
    explicit DumpReader(const std::string &path) : input_(path), path_(path) {
        if (!input_) throw std::runtime_error("cannot open trajectory: " + path);
    }

    bool next(DumpFrame &frame, long long expected_atoms) {
        std::string line;
        while (std::getline(input_, line) && trim(line).empty()) {}
        if (!input_) return false;
        if (trim(line) != "ITEM: TIMESTEP")
            throw std::runtime_error("expected ITEM: TIMESTEP in " + path_);
        if (!std::getline(input_, line)) throw std::runtime_error("missing timestep");
        frame.timestep = std::stoll(trim(line));
        require("ITEM: NUMBER OF ATOMS");
        if (!std::getline(input_, line)) throw std::runtime_error("missing atom count");
        const long long atom_count = std::stoll(trim(line));
        if (atom_count != expected_atoms)
            throw std::runtime_error("trajectory atom count differs from topology data");
        if (!std::getline(input_, line) || !begins_with(trim(line), "ITEM: BOX BOUNDS"))
            throw std::runtime_error("missing BOX BOUNDS");
        frame.box = {};
        read_bounds(frame.box.xlo, frame.box.xhi);
        read_bounds(frame.box.ylo, frame.box.yhi);
        read_bounds(frame.box.zlo, frame.box.zhi);
        frame.box.have_x = frame.box.have_y = frame.box.have_z = true;
        if (!std::getline(input_, line) || !begins_with(trim(line), "ITEM: ATOMS"))
            throw std::runtime_error("missing ATOMS header");
        std::vector<std::string> columns = words(line);
        columns.erase(columns.begin(), columns.begin() + 2);
        const auto column = [&](const std::string &name) {
            const auto found = std::find(columns.begin(), columns.end(), name);
            return found == columns.end() ? -1 :
                static_cast<int>(found - columns.begin());
        };
        const int id_column = column("id");
        const int xu_column = column("xu"), yu_column = column("yu"),
                  zu_column = column("zu");
        const int x_column = column("x"), y_column = column("y"),
                  z_column = column("z");
        const int ix_column = column("ix"), iy_column = column("iy"),
                  iz_column = column("iz");
        const bool have_unwrapped =
            xu_column >= 0 && yu_column >= 0 && zu_column >= 0;
        const bool have_wrapped_images = x_column >= 0 && y_column >= 0 &&
            z_column >= 0 && ix_column >= 0 && iy_column >= 0 && iz_column >= 0;
        if (id_column < 0 || (!have_unwrapped && !have_wrapped_images))
            throw std::runtime_error(
                "trajectory needs id plus xu/yu/zu or x/y/z/ix/iy/iz columns");
        frame.unwrapped.assign(static_cast<std::size_t>(atom_count + 1), {
            std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0});
        for (long long row = 0; row < atom_count; ++row) {
            if (!std::getline(input_, line)) throw std::runtime_error("truncated ATOMS block");
            const std::vector<std::string> values = words(line);
            if (values.size() < columns.size()) throw std::runtime_error("short atom row");
            const long long id = std::stoll(values[static_cast<std::size_t>(id_column)]);
            if (id < 1 || id > atom_count ||
                std::isfinite(frame.unwrapped[static_cast<std::size_t>(id)].x))
                throw std::runtime_error("invalid or duplicate trajectory atom ID");
            Vec3 position;
            if (have_unwrapped) {
                position = {
                    std::stod(values[static_cast<std::size_t>(xu_column)]),
                    std::stod(values[static_cast<std::size_t>(yu_column)]),
                    std::stod(values[static_cast<std::size_t>(zu_column)])};
            } else {
                position = {
                    std::stod(values[static_cast<std::size_t>(x_column)]) +
                        std::stoll(values[static_cast<std::size_t>(ix_column)]) * frame.box.lx(),
                    std::stod(values[static_cast<std::size_t>(y_column)]) +
                        std::stoll(values[static_cast<std::size_t>(iy_column)]) * frame.box.ly(),
                    std::stod(values[static_cast<std::size_t>(z_column)]) +
                        std::stoll(values[static_cast<std::size_t>(iz_column)]) * frame.box.lz()};
            }
            frame.unwrapped[static_cast<std::size_t>(id)] = position;
        }
        return true;
    }

  private:
    void require(const std::string &expected) {
        std::string line;
        if (!std::getline(input_, line) || trim(line) != expected)
            throw std::runtime_error("expected " + expected + " in " + path_);
    }
    void read_bounds(double &lo, double &hi) {
        std::string line;
        if (!std::getline(input_, line)) throw std::runtime_error("truncated BOX BOUNDS");
        std::istringstream fields(line);
        double tilt = 0.0;
        if (!(fields >> lo >> hi)) throw std::runtime_error("invalid box bounds");
        if (fields >> tilt)
            throw std::runtime_error("triclinic trajectories are not currently supported");
    }
    std::ifstream input_;
    std::string path_;
};

int bin_index(double z, const Box &box, int bins, bool periodic_z) {
    if (periodic_z) z = wrap_position({0.0, 0.0, z}, box, true).z;
    int index = static_cast<int>(std::floor((z - box.zlo) / box.lz() * bins));
    if (index < 0) index = 0;
    if (index >= bins) index = bins - 1;
    return index;
}

Vec3 path_midpoint(
    const EffectiveStrand &strand, const DataFile &data,
    const ModelInfo &info) {
    const auto positions = unwrapped_path(strand.atoms, data, info);
    const double half = 0.5 * strand.contour_length;
    double traversed = 0.0;
    for (std::size_t i = 1; i < positions.size(); ++i) {
        const Vec3 segment = positions[i] - positions[i - 1];
        const double length = norm(segment);
        if (traversed + length >= half && length > 0.0) {
            const Vec3 point = positions[i - 1] +
                ((half - traversed) / length) * segment;
            return wrap_position(point, data.box, info.periodic_z());
        }
        traversed += length;
    }
    return wrap_position(positions.back(), data.box, info.periodic_z());
}

Vec3 molecule_center(
    long long molecule, const DataFile &data, const ModelInfo &info) {
    const auto &atoms = data.molecule_atoms[static_cast<std::size_t>(molecule)];
    const Vec3 reference = data.atoms[static_cast<std::size_t>(atoms.front())].position;
    Vec3 sum = reference;
    for (std::size_t i = 1; i < atoms.size(); ++i) {
        const Vec3 position = data.atoms[static_cast<std::size_t>(atoms[i])].position;
        sum += reference + minimum_image_vector(
            position - reference, data.box, info.periodic_z());
    }
    return wrap_position((1.0 / atoms.size()) * sum, data.box, info.periodic_z());
}

std::vector<BinCounts> static_counts(
    const DataFile &data, const ModelInfo &info,
    const ReducedNetwork &network, int bins) {
    std::vector<BinCounts> counts(static_cast<std::size_t>(bins));
    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        const Atom &first = data.atoms[static_cast<std::size_t>(bond.first)];
        const Atom &second = data.atoms[static_cast<std::size_t>(bond.second)];
        const Vec3 delta = minimum_image_vector(
            second.position - first.position, data.box, info.periodic_z());
        const Vec3 midpoint = wrap_position(
            first.position + 0.5 * delta, data.box, info.periodic_z());
        BinCounts &bin = counts[static_cast<std::size_t>(
            bin_index(midpoint.z, data.box, bins, info.periodic_z()))];
        ++bin.crosslink_bonds;
        const int first_component = component_for_molecule(first.molecule, info);
        const int second_component = component_for_molecule(second.molecule, info);
        if (first_component == kStrand || second_component == kStrand)
            ++bin.strand_crosslink_bonds;
        if (first_component == kModerator || second_component == kModerator)
            ++bin.moderator_crosslink_bonds;
    }
    for (std::size_t node = 1; node < network.nodes.size(); ++node)
        ++counts[static_cast<std::size_t>(bin_index(
            network.nodes[node].position.z, data.box, bins,
            info.periodic_z()))].junctions;
    for (const EffectiveStrand &strand : network.strands) {
        Vec3 marker = path_midpoint(strand, data, info);
        BinCounts &bin = counts[static_cast<std::size_t>(bin_index(
            marker.z, data.box, bins, info.periodic_z()))];
        if (strand.status == "active") ++bin.active_strands;
        else if (strand.status == "dangling_loop") ++bin.dangling_loops;
        else if (strand.status == "self_loop") ++bin.self_loops;
        else if (strand.status == "dangling") {
            const long long free_atom = strand.first_node == 0
                ? strand.first_atom : strand.second_atom;
            const double z = data.atoms[static_cast<std::size_t>(free_atom)].position.z;
            ++counts[static_cast<std::size_t>(bin_index(
                z, data.box, bins, info.periodic_z()))].dangling_ends;
        }
    }
    for (const ParentRecord &parent : network.parents) {
        if (parent.state != "isolated" && parent.state != "isolated_ring" &&
            parent.state != "isolated_star" && parent.state != "isolated_grafted") continue;
        const Vec3 center = molecule_center(parent.molecule, data, info);
        ++counts[static_cast<std::size_t>(bin_index(
            center.z, data.box, bins, info.periodic_z()))].isolated_parents;
    }
    return counts;
}

void write_profile(
    const std::filesystem::path &path, const std::vector<BinCounts> &counts,
    const DataFile &data) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "bin\tzlo_A\tzhi_A\tz_center_A\tvolume_A3\tcrosslink_bonds"
        << "\tstrand_crosslink_bonds\tmoderator_crosslink_bonds\tjunctions"
        << "\tactive_strands\tdangling_ends\tdangling_loops\tself_loops"
        << "\tisolated_parents\tcrosslink_density_A-3\tdefect_density_A-3\n";
    const double width = data.box.lz() / counts.size();
    const double volume = data.box.lx() * data.box.ly() * width;
    out << std::setprecision(12);
    for (std::size_t bin = 0; bin < counts.size(); ++bin) {
        const BinCounts &entry = counts[bin];
        const double zlo = data.box.zlo + bin * width;
        const long long defects = entry.dangling_ends + entry.dangling_loops +
            entry.self_loops + entry.isolated_parents;
        out << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
            << zlo + 0.5 * width << '\t' << volume << '\t'
            << entry.crosslink_bonds << '\t' << entry.strand_crosslink_bonds
            << '\t' << entry.moderator_crosslink_bonds << '\t' << entry.junctions
            << '\t' << entry.active_strands << '\t' << entry.dangling_ends
            << '\t' << entry.dangling_loops << '\t' << entry.self_loops
            << '\t' << entry.isolated_parents << '\t'
            << entry.crosslink_bonds / volume << '\t' << defects / volume << '\n';
    }
}

void write_layer_dynamics(
    const std::filesystem::path &path, const Options &options,
    const ModelInfo &info, const DataFile &data, int bins) {
    DumpReader reader(options.trajectory_file);
    DumpFrame origin;
    if (!reader.next(origin, data.declared_atoms))
        throw std::runtime_error("trajectory contains no frames");
    std::vector<int> origin_bin(origin.unwrapped.size(), -1);
    std::vector<long long> bin_counts(static_cast<std::size_t>(bins), 0);
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        if (component_for_molecule(atom.molecule, info) != kStrand) continue;
        const int bin = bin_index(origin.unwrapped[static_cast<std::size_t>(id)].z,
                                  origin.box, bins, info.periodic_z());
        origin_bin[static_cast<std::size_t>(id)] = bin;
        ++bin_counts[static_cast<std::size_t>(bin)];
    }
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "frame\ttimestep\ttime_ns\tbin\tzlo_origin_A\tzhi_origin_A"
        << "\tstrand_beads\tmsd_x_A2\tmsd_y_A2\tmsd_z_A2"
        << "\tmsd_parallel_A2\tmsd_total_A2\tdrift_x_A\tdrift_y_A\tdrift_z_A\n";
    const double width = origin.box.lz() / bins;
    long long frame_index = 0;
    DumpFrame current = origin;
    while (true) {
        if (frame_index % options.frame_stride == 0) {
            Vec3 drift;
            for (long long id = 1; id <= data.declared_atoms; ++id)
                drift += current.unwrapped[static_cast<std::size_t>(id)] -
                         origin.unwrapped[static_cast<std::size_t>(id)];
            drift = (1.0 / data.declared_atoms) * drift;
            std::vector<Vec3> sum(static_cast<std::size_t>(bins));
            for (long long id = 1; id <= data.declared_atoms; ++id) {
                const int bin = origin_bin[static_cast<std::size_t>(id)];
                if (bin < 0) continue;
                const Vec3 displacement =
                    current.unwrapped[static_cast<std::size_t>(id)] -
                    origin.unwrapped[static_cast<std::size_t>(id)] - drift;
                sum[static_cast<std::size_t>(bin)].x += displacement.x * displacement.x;
                sum[static_cast<std::size_t>(bin)].y += displacement.y * displacement.y;
                sum[static_cast<std::size_t>(bin)].z += displacement.z * displacement.z;
            }
            const double time_ns =
                (current.timestep - origin.timestep) * info.timestep_fs * 1.0e-6;
            for (int bin = 0; bin < bins; ++bin) {
                const long long count = bin_counts[static_cast<std::size_t>(bin)];
                const Vec3 msd = count > 0
                    ? (1.0 / count) * sum[static_cast<std::size_t>(bin)] :
                      Vec3{std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()};
                const double zlo = origin.box.zlo + bin * width;
                out << frame_index << '\t' << current.timestep << '\t'
                    << std::setprecision(12) << time_ns << '\t' << bin + 1 << '\t'
                    << zlo << '\t' << zlo + width << '\t' << count << '\t'
                    << msd.x << '\t' << msd.y << '\t' << msd.z << '\t'
                    << msd.x + msd.y << '\t' << msd.x + msd.y + msd.z << '\t'
                    << drift.x << '\t' << drift.y << '\t' << drift.z << '\n';
            }
        }
        ++frame_index;
        if (!reader.next(current, data.declared_atoms)) break;
    }
}

void write_report(
    const std::filesystem::path &path, const Options &options,
    const ModelInfo &info, const DataFile &data, int bins) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "PDMS network distribution and dynamics report\n"
        << "case: " << info.case_name << "\n"
        << "geometry: " << info.geometry << "\n"
        << "z bins: " << bins << "\n"
        << "realized bin width: " << data.box.lz() / bins << " A\n"
        << "static markers: reaction-bond midpoints, reduced-strand midpoints,"
        << " free dangling ends, loop markers, and parent centers\n"
        << "densities use the full lateral box area\n";
    if (options.trajectory_file.empty()) {
        out << "layer dynamics: not requested\n";
    } else {
        out << "trajectory: " << options.trajectory_file << "\n"
            << "layer dynamics selection: all component-1 beads, grouped by first-frame z\n"
            << "MSD convention: displacement from first frame with whole-system COM drift removed\n"
            << "frame stride: " << options.frame_stride << "\n";
    }
    if (info.geometry == "film")
        out << "film guidance: compare x/y or parallel MSD; z is confined\n";
}

void print_help(const char *program) {
    std::cout
        << "Usage: " << program << " <case>.npt_eq <case>.info [options]\n\n"
        << "Options:\n"
        << "  --trajectory FILE    optional dump.msd.lammpstrj for layer MSD\n"
        << "  --bin-width X        target z-bin width in A (default 5)\n"
        << "  --frame-stride N     analyze every Nth trajectory frame (default 1)\n"
        << "  --output-dir PATH    output directory (default analysis_<case>)\n"
        << "  --help               show this help\n";
}

Options parse_options(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        print_help(argv[0]); std::exit(0);
    }
    if (argc < 3) throw std::runtime_error("expected data and info files");
    Options options;
    options.data_file = argv[1];
    options.info_file = argv[2];
    for (int index = 3; index < argc; ++index) {
        const std::string option = argv[index];
        const auto value = [&]() {
            if (++index >= argc) throw std::runtime_error("missing value for " + option);
            return std::string(argv[index]);
        };
        if (option == "--trajectory") options.trajectory_file = value();
        else if (option == "--bin-width") options.bin_width = std::stod(value());
        else if (option == "--frame-stride") options.frame_stride = std::stoll(value());
        else if (option == "--output-dir") options.output_directory = value();
        else if (option == "--help") { print_help(argv[0]); std::exit(0); }
        else throw std::runtime_error("unknown option: " + option);
    }
    if (options.bin_width <= 0.0) throw std::runtime_error("bin width must be positive");
    if (options.frame_stride < 1) throw std::runtime_error("frame stride must be positive");
    return options;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        const DataFile data = parse_data_file(options.data_file, info);
        const ReducedNetwork network = reduce_network(data, info);
        const int bins = std::max(1, static_cast<int>(
            std::ceil(data.box.lz() / options.bin_width)));
        const std::filesystem::path directory = analysis_directory(
            options.data_file, info, options.output_directory);
        pdms_analysis::create_directory(directory);
        const std::string name = safe_case_name(info.case_name);
        write_profile(directory / ("network_z_profile." + name + ".tsv"),
                      static_counts(data, info, network, bins), data);
        if (!options.trajectory_file.empty())
            write_layer_dynamics(
                directory / ("layer_dynamics." + name + ".tsv"),
                options, info, data, bins);
        write_report(directory / ("profile_report." + name + ".txt"),
                     options, info, data, bins);
        std::cout << "Network profiles written to " << directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "network_profile_analyzer: " << error.what() << '\n';
        return 1;
    }
}
