#pragma once

#include "cgkmc_core/config.hpp"

#include <array>
#include <cstdint>
#include <cstddef>
#include <iosfwd>
#include <random>
#include <vector>

namespace cgkmc {

class DumpWriter;

class FenwickTree {
public:
    void reset(std::size_t size);
    void set(std::size_t index, double value);

    [[nodiscard]] double total() const;
    [[nodiscard]] double value(std::size_t index) const;
    [[nodiscard]] std::size_t lower_bound(double prefix_mass) const;

private:
    std::vector<double> tree_;
    std::vector<double> values_;
};

class DynamicIndexSet {
public:
    void reset(std::size_t universe_size);
    void add(std::size_t value);
    void remove(std::size_t value);

    [[nodiscard]] bool contains(std::size_t value) const;
    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] std::size_t at(std::size_t index) const;
    [[nodiscard]] const std::vector<std::size_t>& values() const;

private:
    std::vector<std::size_t> values_;
    std::vector<std::size_t> positions_;
    std::vector<std::uint8_t> present_;
};

class Simulation {
public:
    explicit Simulation(SimulationConfig config, std::uint64_t seed = 0);

    void initialize();
    void step();
    void run(DumpWriter& writer, std::size_t dump_every);
    void run(DumpWriter& writer, std::size_t dump_every, std::ostream& log_output);

    [[nodiscard]] std::size_t step_index() const;
    [[nodiscard]] double time() const;
    [[nodiscard]] double total_energy() const;
    [[nodiscard]] double occupancy() const;
    [[nodiscard]] std::size_t occupied_count() const;
    [[nodiscard]] std::size_t site_count() const;
    [[nodiscard]] std::size_t coordination_number() const;
    [[nodiscard]] double cohesive_energy() const;
    [[nodiscard]] std::size_t solid_interface_count() const;
    [[nodiscard]] std::size_t solvent_interface_count() const;
    [[nodiscard]] std::array<double, 3> bounds() const;
    [[nodiscard]] std::array<double, 3> site_position(std::size_t site_id) const;
    [[nodiscard]] bool is_occupied(std::size_t site_id) const;

    void flip_site_for_testing(std::size_t site_id);
    [[nodiscard]] bool incremental_state_matches_for_testing(double tolerance) const;

private:
    struct PatternNeighbor {
        int dx{};
        int dy{};
        int dz{};
        std::uint32_t basis{};
        double energy{};
    };

    SimulationConfig config_;
    std::mt19937_64 rng_;
    std::size_t site_count_{};
    std::size_t basis_count_{};
    std::size_t coordination_number_{};
    double cohesive_energy_{};
    double surface_density_{};
    double kappa_{};
    std::size_t step_index_{};
    double time_{};
    std::size_t occupied_count_{};
    double total_energy_{};

    std::vector<std::size_t> neighbor_offsets_;
    std::vector<std::uint32_t> neighbor_ids_;
    std::vector<double> neighbor_energies_;
    std::vector<std::uint8_t> occupied_;
    std::vector<int> occupied_neighbor_count_;
    std::vector<double> local_field_;
    DynamicIndexSet solid_interface_;
    DynamicIndexSet solvent_interface_;
    FenwickTree evaporation_weights_;

    void build_neighbor_graph();
    void initialize_state();
    void flip_site(std::size_t site_id);
    void update_interface_membership(std::size_t site_id);
    [[nodiscard]] std::size_t site_id(int ix, int iy, int iz, std::size_t basis) const;
    [[nodiscard]] double interaction_energy(double distance) const;
};

} // namespace cgkmc
