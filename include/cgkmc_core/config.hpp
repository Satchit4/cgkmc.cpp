#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <vector>

namespace cgkmc {

struct LatticeConfig {
    std::array<int, 3> dimensions{};
    std::array<double, 3> lattice_parameters{};
    std::vector<std::array<double, 3>> atomic_basis;

    [[nodiscard]] std::size_t basis_count() const;
    [[nodiscard]] std::size_t site_count() const;
    [[nodiscard]] double density() const;
    [[nodiscard]] double molecular_volume() const;
};

struct InteractionConfig {
    std::vector<double> cutoffs;
    std::vector<double> interaction_energies;
};

struct SolventConfig {
    double beta{};
    double diffusivity{};
    double solubility_limit{};
};

struct GrowthConfig {
    double initial_radius{};
    std::size_t num_steps{};
    std::size_t desired_size{};
};

struct SimulationConfig {
    LatticeConfig lattice;
    InteractionConfig interactions;
    SolventConfig solvent;
    GrowthConfig growth;
};

[[nodiscard]] SimulationConfig load_config(const std::filesystem::path& path);

} // namespace cgkmc
