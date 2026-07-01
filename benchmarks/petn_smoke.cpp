#include "cgkmc_core/config.hpp"
#include "cgkmc_core/dump_writer.hpp"
#include "cgkmc_core/simulation.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>

namespace {

[[nodiscard]] cgkmc::SimulationConfig petn_config(std::size_t steps)
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
    config.growth.num_steps = steps;
    config.growth.desired_size = 40000;
    return config;
}

} // namespace

int main()
{
    cgkmc::Simulation no_dump(petn_config(1000), 0);
    no_dump.initialize();
    for (std::size_t step = 0; step < 1000; ++step) {
        no_dump.step();
    }
    std::cout << "no_dump_steps " << no_dump.step_index() << '\n';
    std::cout << "no_dump_energy " << no_dump.total_energy() << '\n';

    cgkmc::Simulation with_dump(petn_config(1000), 0);
    with_dump.initialize();
    std::ofstream sink("/dev/null");
    cgkmc::DumpWriter writer(sink);
    with_dump.run(writer, 1000);
    std::cout << "dump_every_1000_steps " << with_dump.step_index() << '\n';
    std::cout << "dump_every_1000_energy " << with_dump.total_energy() << '\n';
    return 0;
}
