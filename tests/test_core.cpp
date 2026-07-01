#include "cgkmc_core/config.hpp"
#include "cgkmc_core/dump_writer.hpp"
#include "cgkmc_core/simulation.hpp"

#include <cassert>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

[[nodiscard]] cgkmc::SimulationConfig base_config()
{
    cgkmc::SimulationConfig config;
    config.lattice.dimensions = {5, 5, 5};
    config.lattice.lattice_parameters = {1.0, 1.0, 1.0};
    config.lattice.atomic_basis = {{0.0, 0.0, 0.0}};
    config.interactions.cutoffs = {1.1};
    config.interactions.interaction_energies = {-1.0};
    config.solvent.beta = 38.68;
    config.solvent.diffusivity = 1.0e11;
    config.solvent.solubility_limit = 1.0e-4;
    config.growth.initial_radius = 2.0;
    config.growth.num_steps = 10;
    config.growth.desired_size = 4000;
    return config;
}

void test_coordination_numbers()
{
    {
        cgkmc::SimulationConfig config = base_config();
        config.lattice.atomic_basis = {{0.0, 0.0, 0.0}};
        config.interactions.cutoffs = {1.1};
        cgkmc::Simulation simulation(config);
        simulation.initialize();
        assert(simulation.coordination_number() == 6);
    }
    {
        cgkmc::SimulationConfig config = base_config();
        config.lattice.atomic_basis = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
        config.interactions.cutoffs = {0.9};
        cgkmc::Simulation simulation(config);
        simulation.initialize();
        assert(simulation.coordination_number() == 8);
    }
    {
        cgkmc::SimulationConfig config = base_config();
        config.lattice.atomic_basis = {
            {0.0, 0.0, 0.0},
            {0.5, 0.5, 0.0},
            {0.5, 0.0, 0.5},
            {0.0, 0.5, 0.5},
        };
        config.interactions.cutoffs = {0.8};
        cgkmc::Simulation simulation(config);
        simulation.initialize();
        assert(simulation.coordination_number() == 12);
    }
}

void test_petn_graph()
{
    cgkmc::SimulationConfig config;
    config.lattice.dimensions = {62, 62, 140};
    config.lattice.lattice_parameters = {9.087, 9.087, 6.738};
    config.lattice.atomic_basis = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
    config.interactions.cutoffs = {7.0, 7.5, 9.5};
    config.interactions.interaction_energies = {-0.294, -0.184, -0.002};
    config.solvent.beta = 38.68;
    config.solvent.diffusivity = 1.0e11;
    config.solvent.solubility_limit = 1.0e-4;
    config.growth.initial_radius = 75.0;
    config.growth.num_steps = 10;
    config.growth.desired_size = 40000;

    cgkmc::Simulation simulation(config);
    simulation.initialize();
    assert(simulation.site_count() == 62ULL * 62ULL * 140ULL * 2ULL);
    assert(simulation.coordination_number() == 14);
    assert(std::abs(simulation.cohesive_energy() + 1.034) < 1.0e-9);
}

void test_incremental_updates()
{
    cgkmc::SimulationConfig config = base_config();
    cgkmc::Simulation simulation(config);
    simulation.initialize();
    assert(simulation.incremental_state_matches_for_testing(1.0e-12));

    const std::vector<std::size_t> flips = {0, 31, 62, 63, 64, 12, 124, 1};
    for (const std::size_t site : flips) {
        simulation.flip_site_for_testing(site);
        assert(simulation.incremental_state_matches_for_testing(1.0e-12));
    }
}

void test_dump_format()
{
    cgkmc::SimulationConfig config = base_config();
    config.lattice.dimensions = {3, 3, 3};
    config.growth.initial_radius = 1.1;
    config.growth.num_steps = 1;
    cgkmc::Simulation simulation(config);
    simulation.initialize();

    std::ostringstream output;
    cgkmc::DumpWriter writer(output);
    writer.write_frame(simulation);
    const std::string dump = output.str();
    assert(dump.find("ITEM: TIMESTEP\n0 0.0\n") == 0);
    assert(dump.find("ITEM: NUMBER OF ATOMS\n7\n") != std::string::npos);
    assert(dump.find("ITEM: BOX BOUNDS xy xz xx yy zz\n0.0 3.0 0.0\n") != std::string::npos);
    assert(dump.find("ITEM: ATOMS id type x y z\n") != std::string::npos);
}

void test_small_run()
{
    cgkmc::SimulationConfig config = base_config();
    config.growth.num_steps = 20;
    cgkmc::Simulation simulation(config, 123);
    simulation.initialize();
    for (std::size_t i = 0; i < config.growth.num_steps; ++i) {
        simulation.step();
        assert(simulation.incremental_state_matches_for_testing(1.0e-9));
    }
    assert(simulation.step_index() == config.growth.num_steps);
}

} // namespace

int main()
{
    test_coordination_numbers();
    test_petn_graph();
    test_incremental_updates();
    test_dump_format();
    test_small_run();
    return 0;
}
