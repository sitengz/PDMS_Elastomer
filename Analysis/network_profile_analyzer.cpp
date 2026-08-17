#include "network_common.hpp"

#include <iostream>

namespace {
using namespace pdms_analysis;

struct Options {
    std::string data_file;
    std::string info_file;
    std::string trajectory_file;
    std::string z1_sp_file;
    bool disable_z1 = false;
    std::string output_directory;
    double bin_width = 5.0;
    double wall_fraction = 0.20;
    double core_fraction = 0.20;
    long long frame_stride = 1;
};

constexpr double kAvogadroAngstrom = 0.602214076;

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

struct BinProfile {
    BinCounts topology;
    std::array<long long, kComponentCount> beads{{0, 0, 0, 0}};
    std::array<double, kComponentCount> mass{{0.0, 0.0, 0.0, 0.0}};
    // Strand, crosslinker, and moderator functional sites, respectively.
    std::array<long long, 3> functional_total{{0, 0, 0}};
    std::array<long long, 3> functional_reacted{{0, 0, 0}};

    long long strand_markers = 0;
    long long strand_orientation_count = 0;
    double strand_ree_x2_sum = 0.0;
    double strand_ree_y2_sum = 0.0;
    double strand_ree_z2_sum = 0.0;
    double strand_stretch_sum = 0.0;
    double strand_p2z_sum = 0.0;
    long long active_orientation_count = 0;
    double active_ree_x2_sum = 0.0;
    double active_ree_y2_sum = 0.0;
    double active_ree_z2_sum = 0.0;
    double active_stretch_sum = 0.0;
    double active_p2z_sum = 0.0;

    long long contour_segments = 0;
    double contour_length = 0.0;
    double active_contour_length = 0.0;
    double defect_contour_length = 0.0;
    double segment_p2x_sum = 0.0;
    double segment_p2y_sum = 0.0;
    double segment_p2z_sum = 0.0;
    long long active_contour_segments = 0;
    double active_segment_p2x_sum = 0.0;
    double active_segment_p2y_sum = 0.0;
    double active_segment_p2z_sum = 0.0;

    long long z1_kinks = 0;
    long long z1_segments = 0;
    double z1_primitive_length = 0.0;
    double z1_segment_p2z_sum = 0.0;
};

struct Z1Point {
    Vec3 position;
    double contour_index = 0.0;
    int kink = 0;
};

struct Z1Result {
    Box box;
    std::vector<std::vector<Z1Point>> chains;
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

double mean_or_nan(double sum, long long count) {
    return count > 0 ? sum / static_cast<double>(count) :
        std::numeric_limits<double>::quiet_NaN();
}

double ratio_or_nan(double numerator, double denominator) {
    return denominator != 0.0 ? numerator / denominator :
        std::numeric_limits<double>::quiet_NaN();
}

void add_profile(BinProfile &target, const BinProfile &source) {
    target.topology.crosslink_bonds += source.topology.crosslink_bonds;
    target.topology.strand_crosslink_bonds +=
        source.topology.strand_crosslink_bonds;
    target.topology.moderator_crosslink_bonds +=
        source.topology.moderator_crosslink_bonds;
    target.topology.junctions += source.topology.junctions;
    target.topology.active_strands += source.topology.active_strands;
    target.topology.dangling_ends += source.topology.dangling_ends;
    target.topology.dangling_loops += source.topology.dangling_loops;
    target.topology.self_loops += source.topology.self_loops;
    target.topology.isolated_parents += source.topology.isolated_parents;
    for (int component = 0; component < kComponentCount; ++component) {
        target.beads[static_cast<std::size_t>(component)] +=
            source.beads[static_cast<std::size_t>(component)];
        target.mass[static_cast<std::size_t>(component)] +=
            source.mass[static_cast<std::size_t>(component)];
    }
    for (std::size_t component = 0; component < 3; ++component) {
        target.functional_total[component] += source.functional_total[component];
        target.functional_reacted[component] += source.functional_reacted[component];
    }
#define ADD_FIELD(field) target.field += source.field
    ADD_FIELD(strand_markers);
    ADD_FIELD(strand_orientation_count);
    ADD_FIELD(strand_ree_x2_sum); ADD_FIELD(strand_ree_y2_sum);
    ADD_FIELD(strand_ree_z2_sum); ADD_FIELD(strand_stretch_sum);
    ADD_FIELD(strand_p2z_sum);
    ADD_FIELD(active_orientation_count);
    ADD_FIELD(active_ree_x2_sum); ADD_FIELD(active_ree_y2_sum);
    ADD_FIELD(active_ree_z2_sum); ADD_FIELD(active_stretch_sum);
    ADD_FIELD(active_p2z_sum);
    ADD_FIELD(contour_segments); ADD_FIELD(contour_length);
    ADD_FIELD(active_contour_length); ADD_FIELD(defect_contour_length);
    ADD_FIELD(segment_p2x_sum); ADD_FIELD(segment_p2y_sum);
    ADD_FIELD(segment_p2z_sum); ADD_FIELD(active_contour_segments);
    ADD_FIELD(active_segment_p2x_sum); ADD_FIELD(active_segment_p2y_sum);
    ADD_FIELD(active_segment_p2z_sum);
    ADD_FIELD(z1_kinks); ADD_FIELD(z1_segments);
    ADD_FIELD(z1_primitive_length); ADD_FIELD(z1_segment_p2z_sum);
#undef ADD_FIELD
}

Z1Result read_z1_sp(const std::string &path, const DataFile &data) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Z1+ result: " + path);
    long long chain_count = 0;
    Z1Result result;
    if (!(input >> chain_count) || chain_count < 0 ||
        !(input >> result.box.xhi >> result.box.yhi >> result.box.zhi))
        throw std::runtime_error("invalid Z1+SP header in " + path);
    result.box.have_x = result.box.have_y = result.box.have_z = true;
    const auto close_length = [](double first, double second) {
        return std::fabs(first - second) <=
            1.0e-6 * std::max({1.0, std::fabs(first), std::fabs(second)});
    };
    if (!close_length(result.box.lx(), data.box.lx()) ||
        !close_length(result.box.ly(), data.box.ly()) ||
        !close_length(result.box.lz(), data.box.lz()))
        throw std::runtime_error(
            "Z1+SP box lengths differ from the LAMMPS snapshot; "
            "use unscaled Z1+ coordinates for a physical z profile");
    result.chains.reserve(static_cast<std::size_t>(chain_count));
    for (long long chain = 0; chain < chain_count; ++chain) {
        long long point_count = 0;
        if (!(input >> point_count) || point_count < 2)
            throw std::runtime_error("invalid Z1+ primitive path point count");
        std::vector<Z1Point> points(static_cast<std::size_t>(point_count));
        for (Z1Point &point : points)
            if (!(input >> point.position.x >> point.position.y >> point.position.z >>
                  point.contour_index >> point.kink))
                throw std::runtime_error("truncated Z1+ primitive path data");
        result.chains.push_back(std::move(points));
    }
    std::string trailing;
    if (input >> trailing)
        throw std::runtime_error("unexpected trailing data in Z1+SP result");
    return result;
}

std::set<long long> reacted_functional_atoms(
    const DataFile &data, const ModelInfo &info) {
    std::set<long long> reacted;
    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        reacted.insert(bond.first);
        reacted.insert(bond.second);
    }
    return reacted;
}

std::vector<BinProfile> static_profile(
    const DataFile &data, const ModelInfo &info,
    const ReducedNetwork &network, int bins,
    const Z1Result *z1_result) {
    std::vector<BinProfile> profile(static_cast<std::size_t>(bins));
    const std::set<long long> reacted = reacted_functional_atoms(data, info);

    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        const int component = component_for_molecule(atom.molecule, info);
        BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
            atom.position.z, data.box, bins, info.periodic_z()))];
        ++bin.beads[static_cast<std::size_t>(component)];
        bin.mass[static_cast<std::size_t>(component)] += atom.mass;
    }

    const auto add_declared_sites = [&](int component,
                                        const std::vector<long long> &sites,
                                        std::size_t profile_component) {
        const ComponentInfo &entry = info.components[static_cast<std::size_t>(component)];
        const std::set<long long> declared(sites.begin(), sites.end());
        for (long long molecule = entry.molecule_start;
             molecule > 0 && molecule <= entry.molecule_end; ++molecule) {
            const auto &atoms = data.molecule_atoms[static_cast<std::size_t>(molecule)];
            for (long long atom_id : atoms) {
                if (!declared.count(local_rank(atom_id, atoms))) continue;
                BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
                    data.atoms[static_cast<std::size_t>(atom_id)].position.z,
                    data.box, bins, info.periodic_z()))];
                ++bin.functional_total[profile_component];
                if (reacted.count(atom_id))
                    ++bin.functional_reacted[profile_component];
            }
        }
    };
    add_declared_sites(kStrand, info.reactive_bead_sites, 0);
    add_declared_sites(kCrosslinker, info.crosslinker_reactive_bead_sites, 1);
    const ComponentInfo &moderators = info.components[kModerator];
    for (long long molecule = moderators.molecule_start;
         molecule > 0 && molecule <= moderators.molecule_end; ++molecule) {
        for (long long atom_id :
             data.molecule_atoms[static_cast<std::size_t>(molecule)]) {
            const Atom &atom = data.atoms[static_cast<std::size_t>(atom_id)];
            // The implemented five-bead moderator has a neutral center and
            // four type-2 functional arms.
            if (atom.type != 2) continue;
            BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
                atom.position.z, data.box, bins, info.periodic_z()))];
            ++bin.functional_total[2];
            if (reacted.count(atom_id)) ++bin.functional_reacted[2];
        }
    }

    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        const Atom &first = data.atoms[static_cast<std::size_t>(bond.first)];
        const Atom &second = data.atoms[static_cast<std::size_t>(bond.second)];
        const Vec3 delta = minimum_image_vector(
            second.position - first.position, data.box, info.periodic_z());
        const Vec3 midpoint = wrap_position(
            first.position + 0.5 * delta, data.box, info.periodic_z());
        BinCounts &bin = profile[static_cast<std::size_t>(
            bin_index(midpoint.z, data.box, bins, info.periodic_z()))].topology;
        ++bin.crosslink_bonds;
        const int first_component = component_for_molecule(first.molecule, info);
        const int second_component = component_for_molecule(second.molecule, info);
        if (first_component == kStrand || second_component == kStrand)
            ++bin.strand_crosslink_bonds;
        if (first_component == kModerator || second_component == kModerator)
            ++bin.moderator_crosslink_bonds;
    }
    for (std::size_t node = 1; node < network.nodes.size(); ++node)
        ++profile[static_cast<std::size_t>(bin_index(
            network.nodes[node].position.z, data.box, bins,
            info.periodic_z()))].topology.junctions;
    for (const EffectiveStrand &strand : network.strands) {
        Vec3 marker = path_midpoint(strand, data, info);
        BinProfile &marker_bin = profile[static_cast<std::size_t>(bin_index(
            marker.z, data.box, bins, info.periodic_z()))];
        BinCounts &bin = marker_bin.topology;
        ++marker_bin.strand_markers;
        if (strand.end_to_end_length > 0.0) {
            const double inverse_ree2 =
                1.0 / (strand.end_to_end_length * strand.end_to_end_length);
            marker_bin.strand_ree_x2_sum += strand.end_to_end.x * strand.end_to_end.x;
            marker_bin.strand_ree_y2_sum += strand.end_to_end.y * strand.end_to_end.y;
            marker_bin.strand_ree_z2_sum += strand.end_to_end.z * strand.end_to_end.z;
            marker_bin.strand_p2z_sum +=
                0.5 * (3.0 * strand.end_to_end.z * strand.end_to_end.z * inverse_ree2 - 1.0);
            marker_bin.strand_stretch_sum +=
                ratio_or_nan(strand.end_to_end_length, strand.contour_length);
            ++marker_bin.strand_orientation_count;
        }
        if (strand.status == "active") {
            ++bin.active_strands;
            if (strand.end_to_end_length > 0.0) {
                const double inverse_ree2 =
                    1.0 / (strand.end_to_end_length * strand.end_to_end_length);
                marker_bin.active_ree_x2_sum += strand.end_to_end.x * strand.end_to_end.x;
                marker_bin.active_ree_y2_sum += strand.end_to_end.y * strand.end_to_end.y;
                marker_bin.active_ree_z2_sum += strand.end_to_end.z * strand.end_to_end.z;
                marker_bin.active_p2z_sum +=
                    0.5 * (3.0 * strand.end_to_end.z * strand.end_to_end.z * inverse_ree2 - 1.0);
                marker_bin.active_stretch_sum +=
                    ratio_or_nan(strand.end_to_end_length, strand.contour_length);
                ++marker_bin.active_orientation_count;
            }
        } else if (strand.status == "dangling_loop") ++bin.dangling_loops;
        else if (strand.status == "self_loop") ++bin.self_loops;
        else if (strand.status == "dangling") {
            const long long free_atom = strand.first_node == 0
                ? strand.first_atom : strand.second_atom;
            const double z = data.atoms[static_cast<std::size_t>(free_atom)].position.z;
            ++profile[static_cast<std::size_t>(bin_index(
                z, data.box, bins, info.periodic_z()))].topology.dangling_ends;
        }

        const std::vector<Vec3> positions = unwrapped_path(strand.atoms, data, info);
        for (std::size_t point = 1; point < positions.size(); ++point) {
            const Vec3 segment = positions[point] - positions[point - 1];
            const double length = norm(segment);
            if (!(length > 0.0)) continue;
            const Vec3 midpoint = wrap_position(
                positions[point - 1] + 0.5 * segment,
                data.box, info.periodic_z());
            BinProfile &segment_bin = profile[static_cast<std::size_t>(bin_index(
                midpoint.z, data.box, bins, info.periodic_z()))];
            const double inverse_length2 = 1.0 / (length * length);
            const double p2x = 0.5 * (3.0 * segment.x * segment.x * inverse_length2 - 1.0);
            const double p2y = 0.5 * (3.0 * segment.y * segment.y * inverse_length2 - 1.0);
            const double p2z = 0.5 * (3.0 * segment.z * segment.z * inverse_length2 - 1.0);
            ++segment_bin.contour_segments;
            segment_bin.contour_length += length;
            segment_bin.segment_p2x_sum += p2x;
            segment_bin.segment_p2y_sum += p2y;
            segment_bin.segment_p2z_sum += p2z;
            if (strand.status == "active") {
                segment_bin.active_contour_length += length;
                ++segment_bin.active_contour_segments;
                segment_bin.active_segment_p2x_sum += p2x;
                segment_bin.active_segment_p2y_sum += p2y;
                segment_bin.active_segment_p2z_sum += p2z;
            } else {
                segment_bin.defect_contour_length += length;
            }
        }
    }
    for (const ParentRecord &parent : network.parents) {
        if (parent.state != "isolated" && parent.state != "isolated_ring" &&
            parent.state != "isolated_star" && parent.state != "isolated_grafted") continue;
        const Vec3 center = molecule_center(parent.molecule, data, info);
        ++profile[static_cast<std::size_t>(bin_index(
            center.z, data.box, bins, info.periodic_z()))].topology.isolated_parents;
    }

    if (z1_result != nullptr) {
        for (const auto &chain : z1_result->chains) {
            for (const Z1Point &point : chain) {
                if (point.kink == 0) continue;
                ++profile[static_cast<std::size_t>(bin_index(
                    point.position.z, data.box, bins,
                    info.periodic_z()))].z1_kinks;
            }
            for (std::size_t point = 1; point < chain.size(); ++point) {
                const Vec3 segment =
                    chain[point].position - chain[point - 1].position;
                const double length = norm(segment);
                if (!(length > 0.0)) continue;
                Vec3 midpoint = chain[point - 1].position + 0.5 * segment;
                midpoint = wrap_position(midpoint, data.box, info.periodic_z());
                BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
                    midpoint.z, data.box, bins, info.periodic_z()))];
                ++bin.z1_segments;
                bin.z1_primitive_length += length;
                bin.z1_segment_p2z_sum += 0.5 *
                    (3.0 * segment.z * segment.z / (length * length) - 1.0);
            }
        }
    }
    return profile;
}

void write_profile(
    const std::filesystem::path &path, const std::vector<BinProfile> &profile,
    const DataFile &data, const ModelInfo &info, bool z1_available) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "bin\tzlo_A\tzhi_A\tz_center_A\tvolume_A3\tcrosslink_bonds"
        << "\tstrand_crosslink_bonds\tmoderator_crosslink_bonds\tjunctions"
        << "\tactive_strands\tdangling_ends\tdangling_loops\tself_loops"
        << "\tisolated_parents\tcrosslink_density_A-3\tdefect_density_A-3"
        << "\tmaterial_z_normalized\tdistance_from_nearest_box_wall_A"
        << "\tdistance_from_effective_wall_A"
        << "\tbeads_total\tbeads_strand\tbeads_crosslinker\tbeads_moderator\tbeads_filler"
        << "\tmass_density_total_g_cm-3\tmass_density_strand_g_cm-3"
        << "\tmass_density_crosslinker_g_cm-3\tmass_density_moderator_g_cm-3"
        << "\tmass_density_filler_g_cm-3"
        << "\tfunctional_sites_total\tfunctional_sites_reacted\tlocal_conversion"
        << "\tstrand_sites_total\tstrand_sites_reacted\tstrand_local_conversion"
        << "\tcrosslinker_sites_total\tcrosslinker_sites_reacted"
        << "\tcrosslinker_local_conversion\tmoderator_sites_total"
        << "\tmoderator_sites_reacted\tmoderator_local_conversion"
        << "\tjunction_density_A-3\tactive_strand_density_A-3"
        << "\tdangling_end_density_A-3\tdangling_loop_density_A-3"
        << "\tself_loop_density_A-3\tisolated_parent_density_A-3"
        << "\tcontour_length_A\tactive_contour_length_A\tdefect_contour_length_A"
        << "\tcontour_length_density_A-2\tactive_contour_length_density_A-2"
        << "\tdefect_contour_length_density_A-2\tcontour_segments"
        << "\tsegment_P2_x\tsegment_P2_y\tsegment_P2_z"
        << "\tactive_contour_segments\tactive_segment_P2_x"
        << "\tactive_segment_P2_y\tactive_segment_P2_z"
        << "\tstrand_markers\tstrand_Ree_x2_mean_A2\tstrand_Ree_y2_mean_A2"
        << "\tstrand_Ree_z2_mean_A2\tstrand_Ree_parallel2_mean_A2"
        << "\tstrand_Ree_over_Lc_mean\tstrand_P2_z"
        << "\tactive_Ree_x2_mean_A2\tactive_Ree_y2_mean_A2"
        << "\tactive_Ree_z2_mean_A2\tactive_Ree_parallel2_mean_A2"
        << "\tactive_Ree_over_Lc_mean\tactive_strand_P2_z"
        << "\tz1_data_available\tz1_kinks\tz1_kink_density_A-3"
        << "\tz1_primitive_segments\tz1_primitive_length_A"
        << "\tz1_primitive_length_density_A-2\tz1_primitive_P2_z\n";
    const double width = data.box.lz() / profile.size();
    const double volume = data.box.lx() * data.box.ly() * width;
    const double wall_cutoff = info.geometry == "film"
        ? info.film_wall_cutoff_per_side_angstrom : 0.0;
    double material_thickness = info.geometry == "film"
        ? info.film_thickness_angstrom : data.box.lz();
    if (!(material_thickness > 0.0))
        material_thickness = data.box.lz() - 2.0 * wall_cutoff;
    out << std::setprecision(12);
    for (std::size_t bin = 0; bin < profile.size(); ++bin) {
        const BinProfile &entry = profile[bin];
        const BinCounts &counts = entry.topology;
        const double zlo = data.box.zlo + bin * width;
        const double center = zlo + 0.5 * width;
        const double distance_box = std::min(
            center - data.box.zlo, data.box.zhi - center);
        const double material_z =
            (center - data.box.zlo - wall_cutoff) / material_thickness;
        const long long defects = counts.dangling_ends + counts.dangling_loops +
            counts.self_loops + counts.isolated_parents;
        long long beads_total = 0;
        double mass_total = 0.0;
        long long sites_total = 0, sites_reacted = 0;
        for (int component = 0; component < kComponentCount; ++component) {
            beads_total += entry.beads[static_cast<std::size_t>(component)];
            mass_total += entry.mass[static_cast<std::size_t>(component)];
        }
        for (std::size_t component = 0; component < 3; ++component) {
            sites_total += entry.functional_total[component];
            sites_reacted += entry.functional_reacted[component];
        }
        const auto density = [&](double value) { return value / volume; };
        const auto mass_density = [&](double mass) {
            return mass / (kAvogadroAngstrom * volume);
        };
        const double z1_kink_density = z1_available
            ? density(entry.z1_kinks) : std::numeric_limits<double>::quiet_NaN();
        const double z1_length_density = z1_available
            ? density(entry.z1_primitive_length) :
              std::numeric_limits<double>::quiet_NaN();
        const double z1_p2z = z1_available
            ? mean_or_nan(entry.z1_segment_p2z_sum, entry.z1_segments) :
              std::numeric_limits<double>::quiet_NaN();
        out << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
            << center << '\t' << volume << '\t'
            << counts.crosslink_bonds << '\t' << counts.strand_crosslink_bonds
            << '\t' << counts.moderator_crosslink_bonds << '\t' << counts.junctions
            << '\t' << counts.active_strands << '\t' << counts.dangling_ends
            << '\t' << counts.dangling_loops << '\t' << counts.self_loops
            << '\t' << counts.isolated_parents << '\t'
            << density(counts.crosslink_bonds) << '\t' << density(defects) << '\t'
            << material_z << '\t' << distance_box << '\t'
            << distance_box - wall_cutoff << '\t'
            << beads_total;
        for (long long value : entry.beads) out << '\t' << value;
        out << '\t' << mass_density(mass_total);
        for (double mass : entry.mass) out << '\t' << mass_density(mass);
        out << '\t' << sites_total << '\t' << sites_reacted << '\t'
            << ratio_or_nan(sites_reacted, sites_total);
        for (std::size_t component = 0; component < 3; ++component)
            out << '\t' << entry.functional_total[component] << '\t'
                << entry.functional_reacted[component] << '\t'
                << ratio_or_nan(entry.functional_reacted[component],
                                entry.functional_total[component]);
        out << '\t' << density(counts.junctions)
            << '\t' << density(counts.active_strands)
            << '\t' << density(counts.dangling_ends)
            << '\t' << density(counts.dangling_loops)
            << '\t' << density(counts.self_loops)
            << '\t' << density(counts.isolated_parents)
            << '\t' << entry.contour_length
            << '\t' << entry.active_contour_length
            << '\t' << entry.defect_contour_length
            << '\t' << density(entry.contour_length)
            << '\t' << density(entry.active_contour_length)
            << '\t' << density(entry.defect_contour_length)
            << '\t' << entry.contour_segments
            << '\t' << mean_or_nan(entry.segment_p2x_sum, entry.contour_segments)
            << '\t' << mean_or_nan(entry.segment_p2y_sum, entry.contour_segments)
            << '\t' << mean_or_nan(entry.segment_p2z_sum, entry.contour_segments)
            << '\t' << entry.active_contour_segments
            << '\t' << mean_or_nan(entry.active_segment_p2x_sum,
                                    entry.active_contour_segments)
            << '\t' << mean_or_nan(entry.active_segment_p2y_sum,
                                    entry.active_contour_segments)
            << '\t' << mean_or_nan(entry.active_segment_p2z_sum,
                                    entry.active_contour_segments)
            << '\t' << entry.strand_markers
            << '\t' << mean_or_nan(entry.strand_ree_x2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_ree_y2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_ree_z2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_ree_x2_sum + entry.strand_ree_y2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_stretch_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_p2z_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_x2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_y2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_z2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_x2_sum + entry.active_ree_y2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_stretch_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_p2z_sum,
                                    entry.active_orientation_count)
            << '\t' << (z1_available ? 1 : 0)
            << '\t' << entry.z1_kinks << '\t' << z1_kink_density
            << '\t' << entry.z1_segments << '\t' << entry.z1_primitive_length
            << '\t' << z1_length_density << '\t' << z1_p2z << '\n';
    }
}

void write_folded_profile(
    const std::filesystem::path &path, const std::vector<BinProfile> &profile,
    const DataFile &data, const ModelInfo &info, bool z1_available) {
    const std::size_t folded_bins = (profile.size() + 1) / 2;
    std::vector<BinProfile> folded(folded_bins);
    std::vector<int> multiplicity(folded_bins, 0);
    for (std::size_t bin = 0; bin < profile.size(); ++bin) {
        const std::size_t folded_bin = std::min(bin, profile.size() - 1 - bin);
        add_profile(folded[folded_bin], profile[bin]);
        ++multiplicity[folded_bin];
    }
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "folded_bin\tdistance_box_lo_A\tdistance_box_hi_A"
        << "\tdistance_box_center_A\tdistance_effective_wall_center_A"
        << "\tnormalized_distance_effective_wall\tpaired_full_bins\tvolume_A3"
        << "\tmass_density_total_g_cm-3\tcrosslink_density_A-3"
        << "\tjunction_density_A-3\tactive_strand_density_A-3"
        << "\tdefect_density_A-3\tlocal_conversion"
        << "\tactive_contour_length_density_A-2\tdefect_contour_length_density_A-2"
        << "\tsegment_P2_z\tactive_segment_P2_z\tactive_strand_P2_z"
        << "\tactive_Ree_parallel2_mean_A2\tactive_Ree_z2_mean_A2"
        << "\tz1_data_available\tz1_kinks\tz1_kink_density_A-3"
        << "\tz1_primitive_length_density_A-2\tz1_primitive_P2_z\n";
    const double width = data.box.lz() / profile.size();
    const double area = data.box.lx() * data.box.ly();
    const double wall_cutoff = info.geometry == "film"
        ? info.film_wall_cutoff_per_side_angstrom : 0.0;
    double material_thickness = info.geometry == "film"
        ? info.film_thickness_angstrom : data.box.lz();
    if (!(material_thickness > 0.0))
        material_thickness = data.box.lz() - 2.0 * wall_cutoff;
    out << std::setprecision(12);
    for (std::size_t bin = 0; bin < folded.size(); ++bin) {
        const BinProfile &entry = folded[bin];
        const double distance_lo = bin * width;
        const double distance_hi = std::min((bin + 1) * width, 0.5 * data.box.lz());
        const double center = 0.5 * (distance_lo + distance_hi);
        const double volume = area * width * multiplicity[bin];
        double mass = 0.0;
        long long total_sites = 0, reacted_sites = 0;
        for (double value : entry.mass) mass += value;
        for (std::size_t component = 0; component < 3; ++component) {
            total_sites += entry.functional_total[component];
            reacted_sites += entry.functional_reacted[component];
        }
        const long long defects = entry.topology.dangling_ends +
            entry.topology.dangling_loops + entry.topology.self_loops +
            entry.topology.isolated_parents;
        const auto density = [&](double value) { return value / volume; };
        out << bin + 1 << '\t' << distance_lo << '\t' << distance_hi << '\t'
            << center << '\t' << center - wall_cutoff << '\t'
            << 2.0 * (center - wall_cutoff) / material_thickness << '\t'
            << multiplicity[bin] << '\t' << volume << '\t'
            << mass / (kAvogadroAngstrom * volume) << '\t'
            << density(entry.topology.crosslink_bonds) << '\t'
            << density(entry.topology.junctions) << '\t'
            << density(entry.topology.active_strands) << '\t'
            << density(defects) << '\t'
            << ratio_or_nan(reacted_sites, total_sites) << '\t'
            << density(entry.active_contour_length) << '\t'
            << density(entry.defect_contour_length) << '\t'
            << mean_or_nan(entry.segment_p2z_sum, entry.contour_segments) << '\t'
            << mean_or_nan(entry.active_segment_p2z_sum,
                            entry.active_contour_segments) << '\t'
            << mean_or_nan(entry.active_p2z_sum,
                            entry.active_orientation_count) << '\t'
            << mean_or_nan(entry.active_ree_x2_sum + entry.active_ree_y2_sum,
                            entry.active_orientation_count) << '\t'
            << mean_or_nan(entry.active_ree_z2_sum,
                            entry.active_orientation_count) << '\t'
            << (z1_available ? 1 : 0) << '\t' << entry.z1_kinks << '\t'
            << (z1_available ? density(entry.z1_kinks) :
                std::numeric_limits<double>::quiet_NaN()) << '\t'
            << (z1_available ? density(entry.z1_primitive_length) :
                std::numeric_limits<double>::quiet_NaN()) << '\t'
            << (z1_available ? mean_or_nan(entry.z1_segment_p2z_sum,
                                           entry.z1_segments) :
                std::numeric_limits<double>::quiet_NaN()) << '\n';
    }
}

struct RegionProfile {
    BinProfile profile;
    double volume = 0.0;
    long long bins = 0;
};

void write_profile_summary(
    const std::filesystem::path &path, const std::vector<BinProfile> &profile,
    const DataFile &data, const ModelInfo &info, const Options &options,
    bool z1_available) {
    const double width = data.box.lz() / profile.size();
    const double bin_volume = data.box.lx() * data.box.ly() * width;
    const double wall_cutoff = info.geometry == "film"
        ? info.film_wall_cutoff_per_side_angstrom : 0.0;
    double material_thickness = info.geometry == "film"
        ? info.film_thickness_angstrom : data.box.lz();
    if (!(material_thickness > 0.0))
        material_thickness = data.box.lz() - 2.0 * wall_cutoff;
    RegionProfile lower, upper, core;
    const auto include = [&](RegionProfile &region, const BinProfile &entry) {
        add_profile(region.profile, entry);
        region.volume += bin_volume;
        ++region.bins;
    };
    for (std::size_t bin = 0; bin < profile.size(); ++bin) {
        const double center = data.box.zlo + (bin + 0.5) * width;
        const double normalized =
            (center - data.box.zlo - wall_cutoff) / material_thickness;
        if (normalized >= 0.0 && normalized < options.wall_fraction)
            include(lower, profile[bin]);
        if (normalized <= 1.0 && normalized > 1.0 - options.wall_fraction)
            include(upper, profile[bin]);
        if (std::fabs(normalized - 0.5) <= 0.5 * options.core_fraction)
            include(core, profile[bin]);
    }
    RegionProfile wall = lower;
    add_profile(wall.profile, upper.profile);
    wall.volume += upper.volume;
    wall.bins += upper.bins;

    const auto density = [](const RegionProfile &region, double value) {
        return region.volume > 0.0 ? value / region.volume :
            std::numeric_limits<double>::quiet_NaN();
    };
    const auto total_mass = [](const RegionProfile &region) {
        return std::accumulate(region.profile.mass.begin(),
                               region.profile.mass.end(), 0.0);
    };
    const auto total_sites = [](const RegionProfile &region) {
        return std::accumulate(region.profile.functional_total.begin(),
                               region.profile.functional_total.end(), 0LL);
    };
    const auto reacted_sites = [](const RegionProfile &region) {
        return std::accumulate(region.profile.functional_reacted.begin(),
                               region.profile.functional_reacted.end(), 0LL);
    };
    const auto defect_count = [](const RegionProfile &region) {
        const BinCounts &topology = region.profile.topology;
        return topology.dangling_ends + topology.dangling_loops +
            topology.self_loops + topology.isolated_parents;
    };

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "metric\tlower_wall\tupper_wall\tcombined_wall\tcore"
        << "\twall_over_core\tupper_minus_lower_relative\tunits\n";
    out << std::setprecision(12);
    const auto line = [&](const std::string &name, const std::string &units,
                          const auto &value) {
        const double lo = value(lower), hi = value(upper);
        const double wall_value = value(wall), core_value = value(core);
        const double wall_mean = 0.5 * (lo + hi);
        out << name << '\t' << lo << '\t' << hi << '\t' << wall_value
            << '\t' << core_value << '\t'
            << ratio_or_nan(wall_value, core_value) << '\t'
            << ratio_or_nan(hi - lo, wall_mean) << '\t' << units << '\n';
    };
    line("mass_density", "g/cm^3", [&](const RegionProfile &region) {
        return density(region, total_mass(region)) / kAvogadroAngstrom;
    });
    line("crosslink_bond_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, region.profile.topology.crosslink_bonds);
    });
    line("junction_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, region.profile.topology.junctions);
    });
    line("active_strand_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, region.profile.topology.active_strands);
    });
    line("defect_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, defect_count(region));
    });
    line("local_conversion", "fraction", [&](const RegionProfile &region) {
        return ratio_or_nan(reacted_sites(region), total_sites(region));
    });
    line("active_contour_length_density", "A^-2",
         [&](const RegionProfile &region) {
        return density(region, region.profile.active_contour_length);
    });
    line("defect_contour_length_density", "A^-2",
         [&](const RegionProfile &region) {
        return density(region, region.profile.defect_contour_length);
    });
    line("segment_P2_z", "dimensionless", [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.segment_p2z_sum,
                           region.profile.contour_segments);
    });
    line("active_segment_P2_z", "dimensionless",
         [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_segment_p2z_sum,
                           region.profile.active_contour_segments);
    });
    line("active_strand_P2_z", "dimensionless",
         [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_p2z_sum,
                           region.profile.active_orientation_count);
    });
    line("active_Ree_parallel2_mean", "A^2",
         [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_ree_x2_sum +
                           region.profile.active_ree_y2_sum,
                           region.profile.active_orientation_count);
    });
    line("active_Ree_z2_mean", "A^2", [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_ree_z2_sum,
                           region.profile.active_orientation_count);
    });
    if (z1_available) {
        line("z1_kink_density", "A^-3", [&](const RegionProfile &region) {
            return density(region, region.profile.z1_kinks);
        });
        line("z1_primitive_length_density", "A^-2",
             [&](const RegionProfile &region) {
            return density(region, region.profile.z1_primitive_length);
        });
        line("z1_primitive_P2_z", "dimensionless",
             [&](const RegionProfile &region) {
            return mean_or_nan(region.profile.z1_segment_p2z_sum,
                               region.profile.z1_segments);
        });
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
    const ModelInfo &info, const DataFile &data, int bins,
    const std::string &z1_sp_file) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "PDMS network distribution and dynamics report\n"
        << "case: " << info.case_name << "\n"
        << "geometry: " << info.geometry << "\n"
        << "z bins: " << bins << "\n"
        << "realized bin width: " << data.box.lz() / bins << " A\n"
        << "static markers: reaction-bond midpoints, reduced-strand midpoints,"
        << " free dangling ends, loop markers, and parent centers\n"
        << "densities use the full lateral box area and physical bin volume\n"
        << "local conversion: reacted functional sites divided by declared"
        << " functional sites at their final z positions\n"
        << "crosslink bond position: periodicity-corrected final bond midpoint\n"
        << "contour profile: reduced-strand segment length assigned by segment midpoint\n"
        << "defect contour: every reduced strand not classified active\n"
        << "segment orientation: P2_alpha=(3*u_alpha^2-1)/2\n"
        << "strand conformation: assigned by contour midpoint\n";
    if (info.geometry == "film") {
        double thickness = info.film_thickness_angstrom;
        if (!(thickness > 0.0))
            thickness = data.box.lz() -
                2.0 * info.film_wall_cutoff_per_side_angstrom;
        out << "film box Lz: " << data.box.lz() << " A\n"
            << "film wall cutoff per side: "
            << info.film_wall_cutoff_per_side_angstrom << " A\n"
            << "nominal wall-free material thickness: " << thickness << " A\n"
            << "material coordinate: (z-zlo-wall_cutoff)/material_thickness\n"
            << "summary wall fraction per side: " << options.wall_fraction << "\n"
            << "summary centered core fraction: " << options.core_fraction << "\n";
    } else {
        out << "bulk profile guidance: z is periodic and has no physical wall;"
            << " folded and boundary/core outputs are uniformity controls only\n";
    }
    if (z1_sp_file.empty()) {
        out << "Z1+ profile: unavailable; z1 density/orientation columns are NaN\n";
    } else {
        out << "Z1+ primitive-path result: " << z1_sp_file << "\n"
            << "Z1+ kink definition: primitive-path points with nonzero kink flag\n"
            << "Z1+ scaling: none; result box lengths must match the snapshot\n";
    }
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
        << "  --z1-sp FILE         override the default Z1+SP.dat path\n"
        << "  --no-z1              disable default Z1+SP.dat auto-detection\n"
        << "  --bin-width X        target z-bin width in A (default 5)\n"
        << "  --wall-fraction X    material fraction at each wall for summary (default 0.20)\n"
        << "  --core-fraction X    centered material fraction for summary (default 0.20)\n"
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
        else if (option == "--z1-sp" || option == "--z1-results")
            options.z1_sp_file = value();
        else if (option == "--no-z1") options.disable_z1 = true;
        else if (option == "--bin-width") options.bin_width = std::stod(value());
        else if (option == "--wall-fraction")
            options.wall_fraction = std::stod(value());
        else if (option == "--core-fraction")
            options.core_fraction = std::stod(value());
        else if (option == "--frame-stride") options.frame_stride = std::stoll(value());
        else if (option == "--output-dir") options.output_directory = value();
        else if (option == "--help") { print_help(argv[0]); std::exit(0); }
        else throw std::runtime_error("unknown option: " + option);
    }
    if (options.bin_width <= 0.0) throw std::runtime_error("bin width must be positive");
    if (!(options.wall_fraction > 0.0 && options.wall_fraction < 0.5))
        throw std::runtime_error("wall fraction must be between 0 and 0.5");
    if (!(options.core_fraction > 0.0 && options.core_fraction <= 1.0))
        throw std::runtime_error("core fraction must be between 0 and 1");
    if (options.disable_z1 && !options.z1_sp_file.empty())
        throw std::runtime_error("--no-z1 cannot be combined with --z1-sp");
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
        std::string z1_sp_file = options.z1_sp_file;
        if (z1_sp_file.empty() && !options.disable_z1) {
            const std::filesystem::path candidate = directory / "Z1+SP.dat";
            if (std::filesystem::exists(candidate)) z1_sp_file = candidate.string();
        }
        Z1Result z1_result;
        const Z1Result *z1 = nullptr;
        if (!z1_sp_file.empty()) {
            z1_result = read_z1_sp(z1_sp_file, data);
            z1 = &z1_result;
        }
        const std::vector<BinProfile> profile =
            static_profile(data, info, network, bins, z1);
        write_profile(directory / ("network_z_profile." + name + ".tsv"),
                      profile, data, info, z1 != nullptr);
        write_folded_profile(
            directory / ("network_z_profile_folded." + name + ".tsv"),
            profile, data, info, z1 != nullptr);
        write_profile_summary(
            directory / ("network_z_profile_summary." + name + ".tsv"),
            profile, data, info, options, z1 != nullptr);
        if (!options.trajectory_file.empty())
            write_layer_dynamics(
                directory / ("layer_dynamics." + name + ".tsv"),
                options, info, data, bins);
        write_report(directory / ("profile_report." + name + ".txt"),
                     options, info, data, bins, z1_sp_file);
        std::cout << "Network profiles written to " << directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "network_profile_analyzer: " << error.what() << '\n';
        return 1;
    }
}
