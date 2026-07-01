#include "cgkmc_core/config.hpp"
#include "cgkmc_core/dump_writer.hpp"
#include "cgkmc_core/simulation.hpp"

#include <fstream>
#include <string>

int main(int argc, char** argv)
{
    const cgkmc::SimulationConfig config = cgkmc::load_config(argv[1]);
    const std::size_t dump_every = static_cast<std::size_t>(std::stoull(argv[3]));

    cgkmc::Simulation simulation(config, 0);
    simulation.initialize();

    std::ofstream dump_file(argv[2]);
    cgkmc::DumpWriter writer(dump_file);
    if (argc > 4) {
        std::ofstream log_file(argv[4]);
        simulation.run(writer, dump_every, log_file);
    } else {
        simulation.run(writer, dump_every);
    }
    return 0;
}
