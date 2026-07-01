#include "cgkmc_core/dump_writer.hpp"

#include "cgkmc_core/simulation.hpp"

#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

namespace cgkmc {
namespace {

[[nodiscard]] std::string pythonish_float(double value)
{
    std::ostringstream stream;
    stream << std::setprecision(15) << std::defaultfloat << value;
    std::string text = stream.str();
    if (text.find('.') == std::string::npos && text.find('e') == std::string::npos
        && text.find('E') == std::string::npos) {
        text += ".0";
    }
    return text;
}

} // namespace

DumpWriter::DumpWriter(std::ostream& output)
    : output_(output)
{
}

void DumpWriter::write_frame(const Simulation& simulation)
{
    const auto box = simulation.bounds();
    output_ << "ITEM: TIMESTEP\n";
    output_ << simulation.step_index() << ' ' << pythonish_float(simulation.time()) << '\n';
    output_ << "ITEM: NUMBER OF ATOMS\n";
    output_ << simulation.occupied_count() << '\n';
    output_ << "ITEM: BOX BOUNDS xy xz xx yy zz\n";
    output_ << "0.0 " << pythonish_float(box[0]) << " 0.0\n";
    output_ << "0.0 " << pythonish_float(box[1]) << " 0.0\n";
    output_ << "0.0 " << pythonish_float(box[2]) << " 0.0\n";
    output_ << "ITEM: ATOMS id type x y z\n";

    for (std::size_t site = 0; site < simulation.site_count(); ++site) {
        if (simulation.is_occupied(site)) {
            const auto position = simulation.site_position(site);
            output_ << site << " 1 " << std::fixed << std::setprecision(4)
                    << position[0] << ' ' << position[1] << ' ' << position[2]
                    << std::defaultfloat << '\n';
        }
    }
}

} // namespace cgkmc
