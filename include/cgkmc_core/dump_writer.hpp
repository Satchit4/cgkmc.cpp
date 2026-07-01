#pragma once

#include <iosfwd>

namespace cgkmc {

class Simulation;

class DumpWriter {
public:
    explicit DumpWriter(std::ostream& output);

    void write_frame(const Simulation& simulation);

private:
    std::ostream& output_;
};

} // namespace cgkmc
