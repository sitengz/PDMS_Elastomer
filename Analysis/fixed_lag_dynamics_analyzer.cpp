#include "network_common.hpp"

#include <iostream>

namespace {
using namespace pdms_analysis;

struct Options {
    std::string data_file;
    std::string info_file;
    std::string msd_trajectory_file;
    std::string debye_waller_trajectory_file;
    std::string output_directory;
    double msd_lag_ns = 10.0;
    double debye_waller_lag_ps = 10.0;
    double bin_width = 5.0;
    long long origin_stride = 1;
    bool recenter_film = true;
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
        if (!std::getline(input_, line))
            throw std::runtime_error("missing timestep in " + path_);
        frame.timestep = std::stoll(trim(line));
        require("ITEM: NUMBER OF ATOMS");
        if (!std::getline(input_, line))
            throw std::runtime_error("missing atom count in " + path_);
        const long long atom_count = std::stoll(trim(line));
        if (atom_count != expected_atoms)
            throw std::runtime_error(
                "trajectory atom count differs from topology data");
        if (!std::getline(input_, line) ||
            !begins_with(trim(line), "ITEM: BOX BOUNDS"))
            throw std::runtime_error("missing BOX BOUNDS in " + path_);
        frame.box = {};
        read_bounds(frame.box.xlo, frame.box.xhi);
        read_bounds(frame.box.ylo, frame.box.yhi);
        read_bounds(frame.box.zlo, frame.box.zhi);
        frame.box.have_x = frame.box.have_y = frame.box.have_z = true;
        if (!std::getline(input_, line) ||
            !begins_with(trim(line), "ITEM: ATOMS"))
            throw std::runtime_error("missing ATOMS header in " + path_);
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
            if (!std::getline(input_, line))
                throw std::runtime_error("truncated ATOMS block in " + path_);
            const std::vector<std::string> values = words(line);
            if (values.size() < columns.size())
                throw std::runtime_error("short atom row in " + path_);
            const long long id =
                std::stoll(values[static_cast<std::size_t>(id_column)]);
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
                        std::stoll(values[static_cast<std::size_t>(ix_column)]) *
                            frame.box.lx(),
                    std::stod(values[static_cast<std::size_t>(y_column)]) +
                        std::stoll(values[static_cast<std::size_t>(iy_column)]) *
                            frame.box.ly(),
                    std::stod(values[static_cast<std::size_t>(z_column)]) +
                        std::stoll(values[static_cast<std::size_t>(iz_column)]) *
                            frame.box.lz()};
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
        if (!std::getline(input_, line))
            throw std::runtime_error("truncated BOX BOUNDS in " + path_);
        std::istringstream fields(line);
        double tilt = 0.0;
        if (!(fields >> lo >> hi))
            throw std::runtime_error("invalid box bounds in " + path_);
        if (fields >> tilt)
            throw std::runtime_error(
                "triclinic trajectories are not currently supported");
    }

    std::ifstream input_;
    std::string path_;
};

int bin_index(double z, const Box &box, int bins, bool periodic_z) {
    if (periodic_z) z = wrap_position({0.0, 0.0, z}, box, true).z;
    int index = static_cast<int>(
        std::floor((z - box.zlo) / box.lz() * bins));
    if (index < 0) index = 0;
    if (index >= bins) index = bins - 1;
    return index;
}

Vec3 system_center(const DumpFrame &frame, long long atom_count) {
    Vec3 center;
    for (long long id = 1; id <= atom_count; ++id)
        center += frame.unwrapped[static_cast<std::size_t>(id)];
    return (1.0 / atom_count) * center;
}

struct ScalarStats {
    long long count = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void add(double value) {
        if (!std::isfinite(value)) return;
        ++count;
        const double delta = value - mean;
        mean += delta / count;
        m2 += delta * (value - mean);
    }

    double value_or_nan() const {
        return count > 0 ? mean : std::numeric_limits<double>::quiet_NaN();
    }

    double sample_sd() const {
        return count > 1 ? std::sqrt(m2 / (count - 1)) :
            std::numeric_limits<double>::quiet_NaN();
    }
};

struct ComponentStats {
    ScalarStats x, y, z, xy, total;

    void add(const Vec3 &value) {
        x.add(value.x);
        y.add(value.y);
        z.add(value.z);
        xy.add(value.x + value.y);
        total.add(value.x + value.y + value.z);
    }
};

struct ComponentSums {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double xy = 0.0;
    double total = 0.0;

    void add_displacement(const Vec3 &sum) {
        x += sum.x;
        y += sum.y;
        z += sum.z;
        xy += sum.x + sum.y;
        total += sum.x + sum.y + sum.z;
    }

    void add_einstein(const Vec3 &sum, double lag_ps) {
        x += sum.x / (2.0 * lag_ps);
        y += sum.y / (2.0 * lag_ps);
        z += sum.z / (2.0 * lag_ps);
        xy += (sum.x + sum.y) / (4.0 * lag_ps);
        total += (sum.x + sum.y + sum.z) / (6.0 * lag_ps);
    }
};

struct LayerStats {
    ScalarStats beads;
    ComponentStats displacement;
    ComponentStats einstein;
    std::array<long long, kComponentCount> component_observations{};
    long long occupied_origins = 0;
    double volume_exposure_A3 = 0.0;
    ComponentSums pooled_displacement;
    ComponentSums pooled_einstein;
};

struct FixedLagResult {
    std::string observable;
    double target_lag_ps = 0.0;
    ScalarStats actual_lag_ps;
    ComponentStats displacement;
    ComponentStats einstein;
    std::vector<LayerStats> layers;
    long long origins = 0;
};

Vec3 einstein_components(const Vec3 &msd, double lag_ps) {
    return {msd.x / (2.0 * lag_ps),
            msd.y / (2.0 * lag_ps),
            msd.z / (2.0 * lag_ps)};
}

void write_component_mean_sd(
    std::ostream &out, const ComponentStats &stats) {
    out << stats.x.value_or_nan() << '\t' << stats.x.sample_sd() << '\t'
        << stats.y.value_or_nan() << '\t' << stats.y.sample_sd() << '\t'
        << stats.z.value_or_nan() << '\t' << stats.z.sample_sd() << '\t'
        << stats.xy.value_or_nan() << '\t' << stats.xy.sample_sd() << '\t'
        << stats.total.value_or_nan() << '\t' << stats.total.sample_sd();
}

void write_component_pooled(
    std::ostream &out, const ComponentSums &sums, long long observations) {
    const double inverse = observations > 0
        ? 1.0 / static_cast<double>(observations)
        : std::numeric_limits<double>::quiet_NaN();
    out << sums.x * inverse << '\t' << sums.y * inverse << '\t'
        << sums.z * inverse << '\t' << sums.xy * inverse << '\t'
        << sums.total * inverse;
}

FixedLagResult analyze_fixed_lag(
    const std::string &trajectory_file, double target_lag_ps,
    bool report_einstein, const std::filesystem::path &origin_path,
    const std::filesystem::path &layer_origin_path,
    const std::filesystem::path &layer_summary_path,
    const std::filesystem::path &layer_pooled_summary_path,
    const Options &options, const ModelInfo &info, const DataFile &data,
    const std::string &observable) {
    const long long target_lag_steps = static_cast<long long>(
        std::llround(target_lag_ps * 1000.0 / info.timestep_fs));
    if (target_lag_steps <= 0)
        throw std::runtime_error("fixed lag rounds to zero timesteps");

    std::vector<long long> strand_ids;
    strand_ids.reserve(static_cast<std::size_t>(data.declared_atoms));
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        if (component_for_molecule(atom.molecule, info) == kStrand)
            strand_ids.push_back(id);
    }
    if (strand_ids.empty())
        throw std::runtime_error("trajectory has no component-1 beads");

    DumpReader origin_reader(trajectory_file);
    DumpReader target_reader(trajectory_file);
    DumpFrame target_after;
    bool have_after = target_reader.next(target_after, data.declared_atoms);
    if (!have_after) throw std::runtime_error("trajectory contains no frames");
    const Box reference_box = target_after.box;
    const int bins = std::max(1, static_cast<int>(
        std::ceil(reference_box.lz() / options.bin_width)));
    const double width = reference_box.lz() / bins;
    const double reference_midplane =
        0.5 * (reference_box.zlo + reference_box.zhi);

    std::ofstream origins_out(origin_path);
    if (!origins_out)
        throw std::runtime_error("cannot write " + origin_path.string());
    std::ofstream layers_out(layer_origin_path);
    if (!layers_out)
        throw std::runtime_error("cannot write " + layer_origin_path.string());

    const std::string value_prefix = report_einstein ? "msd" : "u2";
    origins_out << "origin_index\torigin_frame\torigin_timestep\torigin_time_ps"
        << "\ttarget_timestep\tactual_lag_steps\tactual_lag_ps"
        << "\tstrand_beads\tfilm_recenter_shift_A\t" << value_prefix
        << "_x_A2\t" << value_prefix << "_y_A2\t" << value_prefix
        << "_z_A2\t" << value_prefix << "_xy_A2\t" << value_prefix
        << "_3D_A2";
    if (report_einstein)
        origins_out << "\tD_E_x_A2_per_ps\tD_E_y_A2_per_ps"
            << "\tD_E_z_A2_per_ps\tD_E_xy_A2_per_ps\tD_E_3D_A2_per_ps";
    origins_out << '\n';

    layers_out << "origin_index\torigin_timestep\torigin_time_ps"
        << "\tactual_lag_ps\tbin\tzlo_aligned_A\tzhi_aligned_A"
        << "\tstrand_beads\tfilm_recenter_shift_A\t" << value_prefix
        << "_x_A2\t" << value_prefix << "_y_A2\t" << value_prefix
        << "_z_A2\t" << value_prefix << "_xy_A2\t" << value_prefix
        << "_3D_A2";
    if (report_einstein)
        layers_out << "\tD_E_x_A2_per_ps\tD_E_y_A2_per_ps"
            << "\tD_E_z_A2_per_ps\tD_E_xy_A2_per_ps\tD_E_3D_A2_per_ps";
    layers_out << "\tall_beads\tcrosslinker_beads\tmoderator_beads"
        << "\tfiller_beads\tlayer_volume_A3\tall_bead_density_A-3"
        << "\tstrand_bead_density_A-3\tcrosslinker_bead_density_A-3"
        << "\tmoderator_bead_density_A-3\tfiller_bead_density_A-3";
    layers_out << '\n';
    origins_out << std::setprecision(12);
    layers_out << std::setprecision(12);

    FixedLagResult result;
    result.observable = observable;
    result.target_lag_ps = target_lag_ps;
    result.layers.resize(static_cast<std::size_t>(bins));
    DumpFrame target_before;
    bool have_before = false;
    DumpFrame origin;
    long long origin_frame = -1;
    long long first_timestep = 0;
    while (origin_reader.next(origin, data.declared_atoms)) {
        ++origin_frame;
        if (origin_frame == 0) first_timestep = origin.timestep;
        if (origin_frame % options.origin_stride != 0) continue;
        const long long requested_timestep = origin.timestep + target_lag_steps;
        while (have_after && target_after.timestep < requested_timestep) {
            target_before = std::move(target_after);
            have_before = true;
            have_after = target_reader.next(target_after, data.declared_atoms);
        }
        if (!have_after) break;

        const DumpFrame *target = &target_after;
        if (have_before && target_before.timestep > origin.timestep &&
            requested_timestep - target_before.timestep <=
                target_after.timestep - requested_timestep)
            target = &target_before;
        const long long actual_lag_steps =
            target->timestep - origin.timestep;
        if (actual_lag_steps <= 0) continue;
        const double actual_lag_ps =
            actual_lag_steps * info.timestep_fs * 1.0e-3;

        const Vec3 drift = system_center(*target, data.declared_atoms) -
            system_center(origin, data.declared_atoms);
        double strand_center_z = 0.0;
        for (const long long id : strand_ids)
            strand_center_z +=
                origin.unwrapped[static_cast<std::size_t>(id)].z;
        strand_center_z /= strand_ids.size();
        const double recenter_shift =
            info.geometry == "film" && options.recenter_film
            ? reference_midplane - strand_center_z : 0.0;

        std::vector<Vec3> squared_sum(static_cast<std::size_t>(bins));
        std::vector<long long> bin_counts(static_cast<std::size_t>(bins), 0);
        std::vector<std::array<long long, kComponentCount>> component_counts(
            static_cast<std::size_t>(bins));
        for (long long id = 1; id <= data.declared_atoms; ++id) {
            const Vec3 position =
                origin.unwrapped[static_cast<std::size_t>(id)];
            const int bin = bin_index(
                position.z + recenter_shift, reference_box, bins,
                info.periodic_z());
            const int component = component_for_molecule(
                data.atoms[static_cast<std::size_t>(id)].molecule, info);
            ++component_counts[static_cast<std::size_t>(bin)]
                [static_cast<std::size_t>(component)];
        }
        Vec3 global_sum;
        for (const long long id : strand_ids) {
            const Vec3 origin_position =
                origin.unwrapped[static_cast<std::size_t>(id)];
            const int bin = bin_index(
                origin_position.z + recenter_shift, reference_box, bins,
                info.periodic_z());
            const Vec3 displacement =
                target->unwrapped[static_cast<std::size_t>(id)] -
                origin_position - drift;
            const Vec3 squared{displacement.x * displacement.x,
                               displacement.y * displacement.y,
                               displacement.z * displacement.z};
            squared_sum[static_cast<std::size_t>(bin)] += squared;
            ++bin_counts[static_cast<std::size_t>(bin)];
            global_sum += squared;
        }
        const Vec3 global_msd =
            (1.0 / strand_ids.size()) * global_sum;
        const Vec3 global_d = einstein_components(global_msd, actual_lag_ps);
        const double origin_time_ps =
            (origin.timestep - first_timestep) * info.timestep_fs * 1.0e-3;

        origins_out << result.origins << '\t' << origin_frame << '\t'
            << origin.timestep << '\t' << origin_time_ps << '\t'
            << target->timestep << '\t' << actual_lag_steps << '\t'
            << actual_lag_ps << '\t' << strand_ids.size() << '\t'
            << recenter_shift << '\t' << global_msd.x << '\t'
            << global_msd.y << '\t' << global_msd.z << '\t'
            << global_msd.x + global_msd.y << '\t'
            << global_msd.x + global_msd.y + global_msd.z;
        if (report_einstein)
            origins_out << '\t' << global_d.x << '\t' << global_d.y << '\t'
                << global_d.z << '\t'
                << (global_msd.x + global_msd.y) / (4.0 * actual_lag_ps)
                << '\t'
                << (global_msd.x + global_msd.y + global_msd.z) /
                    (6.0 * actual_lag_ps);
        origins_out << '\n';

        result.actual_lag_ps.add(actual_lag_ps);
        result.displacement.add(global_msd);
        if (report_einstein) {
            result.einstein.x.add(global_d.x);
            result.einstein.y.add(global_d.y);
            result.einstein.z.add(global_d.z);
            result.einstein.xy.add(
                (global_msd.x + global_msd.y) / (4.0 * actual_lag_ps));
            result.einstein.total.add(
                (global_msd.x + global_msd.y + global_msd.z) /
                (6.0 * actual_lag_ps));
        }

        for (int bin = 0; bin < bins; ++bin) {
            const long long count = bin_counts[static_cast<std::size_t>(bin)];
            const auto &counts =
                component_counts[static_cast<std::size_t>(bin)];
            const long long all_count = std::accumulate(
                counts.begin(), counts.end(), 0LL);
            const double layer_volume_A3 =
                origin.box.lx() * origin.box.ly() * width;
            const Vec3 value = count > 0
                ? (1.0 / count) * squared_sum[static_cast<std::size_t>(bin)]
                : Vec3{std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN()};
            const Vec3 d = einstein_components(value, actual_lag_ps);
            const double zlo = reference_box.zlo + bin * width;
            layers_out << result.origins << '\t' << origin.timestep << '\t'
                << origin_time_ps << '\t' << actual_lag_ps << '\t'
                << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
                << count << '\t' << recenter_shift << '\t' << value.x << '\t'
                << value.y << '\t' << value.z << '\t' << value.x + value.y
                << '\t' << value.x + value.y + value.z;
            if (report_einstein)
                layers_out << '\t' << d.x << '\t' << d.y << '\t' << d.z
                    << '\t' << (value.x + value.y) / (4.0 * actual_lag_ps)
                    << '\t' << (value.x + value.y + value.z) /
                        (6.0 * actual_lag_ps);
            layers_out << '\t' << all_count << '\t'
                << counts[static_cast<std::size_t>(kCrosslinker)] << '\t'
                << counts[static_cast<std::size_t>(kModerator)] << '\t'
                << counts[static_cast<std::size_t>(kFiller)] << '\t'
                << layer_volume_A3 << '\t'
                << all_count / layer_volume_A3 << '\t'
                << count / layer_volume_A3 << '\t'
                << counts[static_cast<std::size_t>(kCrosslinker)] /
                    layer_volume_A3 << '\t'
                << counts[static_cast<std::size_t>(kModerator)] /
                    layer_volume_A3 << '\t'
                << counts[static_cast<std::size_t>(kFiller)] /
                    layer_volume_A3;
            layers_out << '\n';

            LayerStats &stats = result.layers[static_cast<std::size_t>(bin)];
            stats.beads.add(static_cast<double>(count));
            stats.volume_exposure_A3 += layer_volume_A3;
            for (int component = 0; component < kComponentCount; ++component) {
                const long long component_count =
                    counts[static_cast<std::size_t>(component)];
                stats.component_observations[static_cast<std::size_t>(component)] +=
                    component_count;
            }
            stats.displacement.add(value);
            if (report_einstein && count > 0) {
                stats.einstein.x.add(d.x);
                stats.einstein.y.add(d.y);
                stats.einstein.z.add(d.z);
                stats.einstein.xy.add(
                    (value.x + value.y) / (4.0 * actual_lag_ps));
                stats.einstein.total.add(
                    (value.x + value.y + value.z) /
                    (6.0 * actual_lag_ps));
            }
            if (count > 0) {
                ++stats.occupied_origins;
                stats.pooled_displacement.add_displacement(
                    squared_sum[static_cast<std::size_t>(bin)]);
                if (report_einstein)
                    stats.pooled_einstein.add_einstein(
                        squared_sum[static_cast<std::size_t>(bin)],
                        actual_lag_ps);
            }
        }
        ++result.origins;
    }
    if (result.origins < 2)
        throw std::runtime_error(
            "fixed-lag analysis found fewer than two valid origins; "
            "the trajectory must extend beyond the requested lag");

    std::ofstream summary(layer_summary_path);
    if (!summary)
        throw std::runtime_error("cannot write " + layer_summary_path.string());
    summary << "bin\tzlo_aligned_A\tzhi_aligned_A\torigins"
        << "\tmean_strand_beads"
        << "\tmean_x_A2\tsd_x_A2\tmean_y_A2\tsd_y_A2"
        << "\tmean_z_A2\tsd_z_A2\tmean_xy_A2\tsd_xy_A2"
        << "\tmean_3D_A2\tsd_3D_A2";
    if (report_einstein)
        summary << "\tmean_D_E_x_A2_per_ps\tsd_D_E_x_A2_per_ps"
            << "\tmean_D_E_y_A2_per_ps\tsd_D_E_y_A2_per_ps"
            << "\tmean_D_E_z_A2_per_ps\tsd_D_E_z_A2_per_ps"
            << "\tmean_D_E_xy_A2_per_ps\tsd_D_E_xy_A2_per_ps"
            << "\tmean_D_E_3D_A2_per_ps\tsd_D_E_3D_A2_per_ps";
    summary << '\n' << std::setprecision(12);
    for (int bin = 0; bin < bins; ++bin) {
        const LayerStats &stats = result.layers[static_cast<std::size_t>(bin)];
        const double zlo = reference_box.zlo + bin * width;
        summary << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
            << stats.displacement.total.count << '\t'
            << stats.beads.value_or_nan() << '\t';
        write_component_mean_sd(summary, stats.displacement);
        if (report_einstein) {
            summary << '\t';
            write_component_mean_sd(summary, stats.einstein);
        }
        summary << '\n';
    }

    std::ofstream pooled(layer_pooled_summary_path);
    if (!pooled)
        throw std::runtime_error(
            "cannot write " + layer_pooled_summary_path.string());
    pooled << "bin\tzlo_aligned_A\tzhi_aligned_A\torigins_total"
        << "\torigins_occupied\torigin_coverage\tmean_layer_volume_A3"
        << "\tall_bead_observations\tstrand_bead_observations"
        << "\tcrosslinker_bead_observations\tmoderator_bead_observations"
        << "\tfiller_bead_observations\tmean_all_beads\tmean_strand_beads"
        << "\tmean_crosslinker_beads\tmean_moderator_beads"
        << "\tmean_filler_beads\tall_bead_density_A-3"
        << "\tstrand_bead_density_A-3\tcrosslinker_bead_density_A-3"
        << "\tmoderator_bead_density_A-3\tfiller_bead_density_A-3"
        << "\tpooled_mean_x_A2\tpooled_mean_y_A2\tpooled_mean_z_A2"
        << "\tpooled_mean_xy_A2\tpooled_mean_3D_A2";
    if (report_einstein)
        pooled << "\tpooled_D_E_x_A2_per_ps\tpooled_D_E_y_A2_per_ps"
            << "\tpooled_D_E_z_A2_per_ps\tpooled_D_E_xy_A2_per_ps"
            << "\tpooled_D_E_3D_A2_per_ps";
    pooled << '\n' << std::setprecision(12);
    for (int bin = 0; bin < bins; ++bin) {
        const LayerStats &stats = result.layers[static_cast<std::size_t>(bin)];
        const double zlo = reference_box.zlo + bin * width;
        const long long strand_observations =
            stats.component_observations[static_cast<std::size_t>(kStrand)];
        const long long all_observations = std::accumulate(
            stats.component_observations.begin(),
            stats.component_observations.end(), 0LL);
        const double inverse_origins = result.origins > 0
            ? 1.0 / static_cast<double>(result.origins)
            : std::numeric_limits<double>::quiet_NaN();
        const double inverse_volume = stats.volume_exposure_A3 > 0.0
            ? 1.0 / stats.volume_exposure_A3
            : std::numeric_limits<double>::quiet_NaN();
        pooled << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
            << result.origins << '\t' << stats.occupied_origins << '\t'
            << stats.occupied_origins * inverse_origins << '\t'
            << stats.volume_exposure_A3 * inverse_origins << '\t'
            << all_observations << '\t' << strand_observations << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kCrosslinker)] << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kModerator)] << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kFiller)] << '\t'
            << all_observations * inverse_origins << '\t'
            << strand_observations * inverse_origins << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kCrosslinker)] * inverse_origins
            << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kModerator)] * inverse_origins
            << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kFiller)] * inverse_origins
            << '\t' << all_observations * inverse_volume << '\t'
            << strand_observations * inverse_volume << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kCrosslinker)] * inverse_volume
            << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kModerator)] * inverse_volume
            << '\t'
            << stats.component_observations[
                   static_cast<std::size_t>(kFiller)] * inverse_volume
            << '\t';
        write_component_pooled(
            pooled, stats.pooled_displacement, strand_observations);
        if (report_einstein) {
            pooled << '\t';
            write_component_pooled(
                pooled, stats.pooled_einstein, strand_observations);
        }
        pooled << '\n';
    }
    return result;
}

void write_global_summary(
    const std::filesystem::path &path, const FixedLagResult &msd,
    const FixedLagResult &dw) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "observable\ttarget_lag_ps\torigins\tmean_actual_lag_ps"
        << "\tsd_actual_lag_ps"
        << "\tmean_x_A2\tsd_x_A2\tmean_y_A2\tsd_y_A2"
        << "\tmean_z_A2\tsd_z_A2\tmean_xy_A2\tsd_xy_A2"
        << "\tmean_3D_A2\tsd_3D_A2"
        << "\tmean_D_E_x_A2_per_ps\tsd_D_E_x_A2_per_ps"
        << "\tmean_D_E_y_A2_per_ps\tsd_D_E_y_A2_per_ps"
        << "\tmean_D_E_z_A2_per_ps\tsd_D_E_z_A2_per_ps"
        << "\tmean_D_E_xy_A2_per_ps\tsd_D_E_xy_A2_per_ps"
        << "\tmean_D_E_3D_A2_per_ps\tsd_D_E_3D_A2_per_ps\n";
    out << std::setprecision(12);
    const auto row = [&](const FixedLagResult &result, bool have_d) {
        out << result.observable << '\t' << result.target_lag_ps << '\t'
            << result.origins << '\t' << result.actual_lag_ps.value_or_nan()
            << '\t' << result.actual_lag_ps.sample_sd() << '\t';
        write_component_mean_sd(out, result.displacement);
        out << '\t';
        if (have_d) {
            write_component_mean_sd(out, result.einstein);
        } else {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            out << nan;
            for (int column = 1; column < 10; ++column) out << '\t' << nan;
        }
        out << '\n';
    };
    row(msd, true);
    row(dw, false);
}

void write_report(
    const std::filesystem::path &path, const Options &options,
    const ModelInfo &info, const FixedLagResult &msd,
    const FixedLagResult &dw) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "PDMS fixed-lag origin-resolved dynamics report\n"
        << "case: " << info.case_name << "\n"
        << "geometry: " << info.geometry << "\n"
        << "selection: all component-1 beads\n"
        << "drift correction: whole-system number-weighted mean displacement\n"
        << "MSD trajectory: " << options.msd_trajectory_file << "\n"
        << "MSD target lag: " << options.msd_lag_ns << " ns\n"
        << "MSD valid origins: " << msd.origins << "\n"
        << "Einstein definitions: D_x=MSD_x/(2 lag), D_xy=MSD_xy/(4 lag),"
        << " D_3D=MSD_3D/(6 lag)\n"
        << "Debye-Waller trajectory: "
        << options.debye_waller_trajectory_file << "\n"
        << "Debye-Waller target lag: " << options.debye_waller_lag_ps
        << " ps\n"
        << "Debye-Waller valid origins: " << dw.origins << "\n"
        << "origin stride: every " << options.origin_stride << " dump frame(s)\n"
        << "target matching: nearest recorded frame, requiring a bracketing"
        << " frame at or after the requested lag\n"
        << "layer assignment: component-1 z position at each time origin\n"
        << "layer pooled means: bead-origin weighted; sparse origins do not"
        << " receive equal weight with dense origins\n"
        << "layer density: all components counted at the same post-release"
        << " time origins used by each fixed-lag observable\n"
        << "film recentering: "
        << (info.geometry == "film" && options.recenter_film ?
            "origin strand midplane aligned to the reference box midplane" :
            "disabled") << "\n"
        << "origin SD warning: origins may be correlated; the reported SD is"
        << " descriptive and raw origin tables should be used for block"
        << " bootstrap or jackknife uncertainty\n"
        << "film guidance: use xy as the primary long-time observable; z and"
        << " 3D values remain available but confined z motion is not ordinary"
        << " diffusion\n";
}

void print_help(const char *program) {
    std::cout
        << "Usage: " << program
        << " <case>.npt_eq <case>.info --msd-trajectory FILE"
        << " --dw-trajectory FILE [options]\n\n"
        << "Options:\n"
        << "  --msd-trajectory FILE long-time layer-dynamics dump (required)\n"
        << "  --dw-trajectory FILE  high-frequency Debye-Waller dump (required)\n"
        << "  --msd-lag-ns X        fixed MSD/Einstein lag (default 10)\n"
        << "  --dw-lag-ps X         fixed u^2 lag (default 10)\n"
        << "  --origin-stride N     use every Nth possible origin (default 1)\n"
        << "  --bin-width X         target aligned z-bin width in A (default 5)\n"
        << "  --no-film-recenter    keep absolute dump z for film origin layers\n"
        << "  --output-dir PATH     output directory (default analysis_<case>)\n"
        << "  --help                show this help\n";
}

Options parse_options(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        print_help(argv[0]);
        std::exit(0);
    }
    if (argc < 3) throw std::runtime_error("expected data and info files");
    Options options;
    options.data_file = argv[1];
    options.info_file = argv[2];
    for (int index = 3; index < argc; ++index) {
        const std::string option = argv[index];
        const auto value = [&]() {
            if (++index >= argc)
                throw std::runtime_error("missing value for " + option);
            return std::string(argv[index]);
        };
        if (option == "--msd-trajectory")
            options.msd_trajectory_file = value();
        else if (option == "--dw-trajectory")
            options.debye_waller_trajectory_file = value();
        else if (option == "--msd-lag-ns")
            options.msd_lag_ns = std::stod(value());
        else if (option == "--dw-lag-ps")
            options.debye_waller_lag_ps = std::stod(value());
        else if (option == "--origin-stride")
            options.origin_stride = std::stoll(value());
        else if (option == "--bin-width")
            options.bin_width = std::stod(value());
        else if (option == "--no-film-recenter")
            options.recenter_film = false;
        else if (option == "--output-dir")
            options.output_directory = value();
        else if (option == "--help") {
            print_help(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + option);
        }
    }
    if (options.msd_trajectory_file.empty())
        throw std::runtime_error("--msd-trajectory is required");
    if (options.debye_waller_trajectory_file.empty())
        throw std::runtime_error("--dw-trajectory is required");
    if (!(options.msd_lag_ns > 0.0))
        throw std::runtime_error("MSD lag must be positive");
    if (!(options.debye_waller_lag_ps > 0.0))
        throw std::runtime_error("Debye-Waller lag must be positive");
    if (options.origin_stride < 1)
        throw std::runtime_error("origin stride must be positive");
    if (!(options.bin_width > 0.0))
        throw std::runtime_error("bin width must be positive");
    return options;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        const DataFile data = parse_data_file(options.data_file, info);
        const std::filesystem::path directory = analysis_directory(
            options.data_file, info, options.output_directory);
        pdms_analysis::create_directory(directory);
        const std::string name = safe_case_name(info.case_name);

        const FixedLagResult msd = analyze_fixed_lag(
            options.msd_trajectory_file, options.msd_lag_ns * 1000.0, true,
            directory / ("fixed_lag_msd_origins." + name + ".tsv"),
            directory / ("fixed_lag_msd_layers." + name + ".tsv"),
            directory / ("fixed_lag_msd_layer_summary." + name + ".tsv"),
            directory /
                ("fixed_lag_msd_layer_pooled_summary." + name + ".tsv"),
            options, info, data, "msd_fixed_lag");
        const FixedLagResult dw = analyze_fixed_lag(
            options.debye_waller_trajectory_file,
            options.debye_waller_lag_ps, false,
            directory / ("fixed_lag_u2_origins." + name + ".tsv"),
            directory / ("fixed_lag_u2_layers." + name + ".tsv"),
            directory / ("fixed_lag_u2_layer_summary." + name + ".tsv"),
            directory /
                ("fixed_lag_u2_layer_pooled_summary." + name + ".tsv"),
            options, info, data, "u2_fixed_lag");
        write_global_summary(
            directory / ("fixed_lag_summary." + name + ".tsv"), msd, dw);
        write_report(
            directory / ("fixed_lag_report." + name + ".txt"),
            options, info, msd, dw);
        std::cout << "Fixed-lag dynamics written to "
                  << directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "fixed_lag_dynamics_analyzer: " << error.what() << '\n';
        return 1;
    }
}
