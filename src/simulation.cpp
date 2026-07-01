#include "cgkmc_core/simulation.hpp"

#include "cgkmc_core/dump_writer.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <ostream>
#include <utility>

namespace cgkmc {
namespace {

[[nodiscard]] int wrap_index(int value, int dimension)
{
    value %= dimension;
    if (value < 0) {
        value += dimension;
    }
    return value;
}

[[nodiscard]] std::size_t ceil_offset(double cutoff, double lattice_parameter)
{
    return static_cast<std::size_t>(std::ceil(cutoff / lattice_parameter)) + 1;
}

} // namespace

void FenwickTree::reset(std::size_t size)
{
    tree_.assign(size + 1, 0.0);
    values_.assign(size, 0.0);
}

void FenwickTree::set(std::size_t index, double value)
{
    const double delta = value - values_[index];
    values_[index] = value;
    for (std::size_t i = index + 1; i < tree_.size(); i += i & (~i + 1)) {
        tree_[i] += delta;
    }
}

double FenwickTree::total() const
{
    double sum = 0.0;
    for (std::size_t i = values_.size(); i > 0; i -= i & (~i + 1)) {
        sum += tree_[i];
    }
    return sum;
}

double FenwickTree::value(std::size_t index) const
{
    return values_[index];
}

std::size_t FenwickTree::lower_bound(double prefix_mass) const
{
    std::size_t index = 0;
    std::size_t bit = 1;
    while ((bit << 1) < tree_.size()) {
        bit <<= 1;
    }
    while (bit != 0) {
        const std::size_t next = index + bit;
        if (next < tree_.size() && tree_[next] <= prefix_mass) {
            index = next;
            prefix_mass -= tree_[next];
        }
        bit >>= 1;
    }
    return index;
}

void DynamicIndexSet::reset(std::size_t universe_size)
{
    values_.clear();
    positions_.assign(universe_size, 0);
    present_.assign(universe_size, 0);
}

void DynamicIndexSet::add(std::size_t value)
{
    if (!present_[value]) {
        present_[value] = 1;
        positions_[value] = values_.size();
        values_.push_back(value);
    }
}

void DynamicIndexSet::remove(std::size_t value)
{
    if (present_[value]) {
        const std::size_t position = positions_[value];
        const std::size_t replacement = values_.back();
        values_[position] = replacement;
        positions_[replacement] = position;
        values_.pop_back();
        present_[value] = 0;
    }
}

bool DynamicIndexSet::contains(std::size_t value) const
{
    return present_[value] != 0;
}

std::size_t DynamicIndexSet::size() const
{
    return values_.size();
}

std::size_t DynamicIndexSet::at(std::size_t index) const
{
    return values_[index];
}

const std::vector<std::size_t>& DynamicIndexSet::values() const
{
    return values_;
}

Simulation::Simulation(SimulationConfig config, std::uint64_t seed)
    : config_(std::move(config))
    , rng_(seed)
{
}

void Simulation::initialize()
{
    site_count_ = config_.lattice.site_count();
    basis_count_ = config_.lattice.basis_count();
    build_neighbor_graph();

    occupied_.assign(site_count_, 0);
    occupied_neighbor_count_.assign(site_count_, 0);
    local_field_.assign(site_count_, 0.0);
    solid_interface_.reset(site_count_);
    solvent_interface_.reset(site_count_);
    evaporation_weights_.reset(site_count_);
    step_index_ = 0;
    time_ = 0.0;

    std::array<double, 3> center{};
    for (std::size_t basis = 0; basis < basis_count_; ++basis) {
        center[0] += config_.lattice.atomic_basis[basis][0];
        center[1] += config_.lattice.atomic_basis[basis][1];
        center[2] += config_.lattice.atomic_basis[basis][2];
    }
    center[0] = ((static_cast<double>(config_.lattice.dimensions[0] - 1) * 0.5) + center[0] / basis_count_)
        * config_.lattice.lattice_parameters[0];
    center[1] = ((static_cast<double>(config_.lattice.dimensions[1] - 1) * 0.5) + center[1] / basis_count_)
        * config_.lattice.lattice_parameters[1];
    center[2] = ((static_cast<double>(config_.lattice.dimensions[2] - 1) * 0.5) + center[2] / basis_count_)
        * config_.lattice.lattice_parameters[2];

    occupied_count_ = 0;
    const double initial_radius_sq = config_.growth.initial_radius * config_.growth.initial_radius;
    for (std::size_t site = 0; site < site_count_; ++site) {
        const auto position = site_position(site);
        const double dx = position[0] - center[0];
        const double dy = position[1] - center[1];
        const double dz = position[2] - center[2];
        if (dx * dx + dy * dy + dz * dz <= initial_radius_sq) {
            occupied_[site] = 1;
            ++occupied_count_;
        }
    }

    initialize_state();
    const double initial_surface_area = 4.0 * std::numbers::pi * config_.growth.initial_radius
        * config_.growth.initial_radius;
    surface_density_ = static_cast<double>(solid_interface_.size()) / initial_surface_area;
    kappa_ = std::cbrt((4.0 / 3.0) * std::numbers::pi * config_.lattice.density())
        * config_.solvent.diffusivity * config_.solvent.solubility_limit
        / (surface_density_ * std::cbrt(static_cast<double>(config_.growth.desired_size)));
}

void Simulation::step()
{
    const double evaporation_total_rate = kappa_ * static_cast<double>(solid_interface_.size());
    const double radius = std::cbrt(
        0.75 * config_.lattice.molecular_volume() * static_cast<double>(occupied_count_) / std::numbers::pi);
    const double adsorption_site_rate = config_.solvent.diffusivity * config_.solvent.solubility_limit
        / (surface_density_ * radius);
    const double adsorption_total_rate = adsorption_site_rate * static_cast<double>(solvent_interface_.size());
    const double total_rate = evaporation_total_rate + adsorption_total_rate;

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double channel_mass = unit(rng_) * total_rate;
    std::size_t event_site = 0;
    if (channel_mass < evaporation_total_rate) {
        event_site = evaporation_weights_.lower_bound(unit(rng_) * evaporation_weights_.total());
    } else {
        std::uniform_int_distribution<std::size_t> solvent_pick(0, solvent_interface_.size() - 1);
        event_site = solvent_interface_.at(solvent_pick(rng_));
    }

    flip_site(event_site);
    std::exponential_distribution<double> residence_time(total_rate);
    time_ += residence_time(rng_);
    ++step_index_;
}

void Simulation::run(DumpWriter& writer, std::size_t dump_every)
{
    while (step_index_ < config_.growth.num_steps) {
        if (step_index_ % dump_every == 0) {
            writer.write_frame(*this);
        }
        step();
    }
}

void Simulation::run(DumpWriter& writer, std::size_t dump_every, std::ostream& log_output)
{
    while (step_index_ < config_.growth.num_steps) {
        if (step_index_ % dump_every == 0) {
            writer.write_frame(*this);
            log_output << "INFO:simulation info TIME=" << time_
                       << " ENERGY=" << total_energy_
                       << " OCCUPANCY=" << occupancy() << '\n';
        }
        step();
    }
}

std::size_t Simulation::step_index() const
{
    return step_index_;
}

double Simulation::time() const
{
    return time_;
}

double Simulation::total_energy() const
{
    return total_energy_;
}

double Simulation::occupancy() const
{
    return static_cast<double>(occupied_count_) / static_cast<double>(site_count_);
}

std::size_t Simulation::occupied_count() const
{
    return occupied_count_;
}

std::size_t Simulation::site_count() const
{
    return site_count_;
}

std::size_t Simulation::coordination_number() const
{
    return coordination_number_;
}

double Simulation::cohesive_energy() const
{
    return cohesive_energy_;
}

std::size_t Simulation::solid_interface_count() const
{
    return solid_interface_.size();
}

std::size_t Simulation::solvent_interface_count() const
{
    return solvent_interface_.size();
}

std::array<double, 3> Simulation::bounds() const
{
    return {
        config_.lattice.lattice_parameters[0] * config_.lattice.dimensions[0],
        config_.lattice.lattice_parameters[1] * config_.lattice.dimensions[1],
        config_.lattice.lattice_parameters[2] * config_.lattice.dimensions[2],
    };
}

std::array<double, 3> Simulation::site_position(std::size_t site_id) const
{
    const std::size_t basis = site_id % basis_count_;
    std::size_t cell = site_id / basis_count_;
    const int iz = static_cast<int>(cell % static_cast<std::size_t>(config_.lattice.dimensions[2]));
    cell /= static_cast<std::size_t>(config_.lattice.dimensions[2]);
    const int iy = static_cast<int>(cell % static_cast<std::size_t>(config_.lattice.dimensions[1]));
    cell /= static_cast<std::size_t>(config_.lattice.dimensions[1]);
    const int ix = static_cast<int>(cell);

    return {
        (static_cast<double>(ix) + config_.lattice.atomic_basis[basis][0]) * config_.lattice.lattice_parameters[0],
        (static_cast<double>(iy) + config_.lattice.atomic_basis[basis][1]) * config_.lattice.lattice_parameters[1],
        (static_cast<double>(iz) + config_.lattice.atomic_basis[basis][2]) * config_.lattice.lattice_parameters[2],
    };
}

bool Simulation::is_occupied(std::size_t site_id) const
{
    return occupied_[site_id] != 0;
}

void Simulation::build_neighbor_graph()
{
    const double max_cutoff = config_.interactions.cutoffs.back();
    const int rx = static_cast<int>(ceil_offset(max_cutoff, config_.lattice.lattice_parameters[0]));
    const int ry = static_cast<int>(ceil_offset(max_cutoff, config_.lattice.lattice_parameters[1]));
    const int rz = static_cast<int>(ceil_offset(max_cutoff, config_.lattice.lattice_parameters[2]));
    const double max_cutoff_sq = max_cutoff * max_cutoff;

    std::vector<std::vector<PatternNeighbor>> pattern(basis_count_);
    for (std::size_t source_basis = 0; source_basis < basis_count_; ++source_basis) {
        for (std::size_t target_basis = 0; target_basis < basis_count_; ++target_basis) {
            for (int dx = -rx; dx <= rx; ++dx) {
                for (int dy = -ry; dy <= ry; ++dy) {
                    for (int dz = -rz; dz <= rz; ++dz) {
                        const double x = (static_cast<double>(dx)
                                             + config_.lattice.atomic_basis[target_basis][0]
                                             - config_.lattice.atomic_basis[source_basis][0])
                            * config_.lattice.lattice_parameters[0];
                        const double y = (static_cast<double>(dy)
                                             + config_.lattice.atomic_basis[target_basis][1]
                                             - config_.lattice.atomic_basis[source_basis][1])
                            * config_.lattice.lattice_parameters[1];
                        const double z = (static_cast<double>(dz)
                                             + config_.lattice.atomic_basis[target_basis][2]
                                             - config_.lattice.atomic_basis[source_basis][2])
                            * config_.lattice.lattice_parameters[2];
                        const double distance_sq = x * x + y * y + z * z;
                        if (distance_sq > 0.0 && distance_sq <= max_cutoff_sq) {
                            pattern[source_basis].push_back(PatternNeighbor{
                                dx,
                                dy,
                                dz,
                                static_cast<std::uint32_t>(target_basis),
                                interaction_energy(std::sqrt(distance_sq)),
                            });
                        }
                    }
                }
            }
        }
    }

    neighbor_offsets_.assign(site_count_ + 1, 0);
    std::size_t edge_count = 0;
    for (int ix = 0; ix < config_.lattice.dimensions[0]; ++ix) {
        for (int iy = 0; iy < config_.lattice.dimensions[1]; ++iy) {
            for (int iz = 0; iz < config_.lattice.dimensions[2]; ++iz) {
                for (std::size_t basis = 0; basis < basis_count_; ++basis) {
                    const std::size_t site = site_id(ix, iy, iz, basis);
                    neighbor_offsets_[site] = edge_count;
                    edge_count += pattern[basis].size();
                }
            }
        }
    }
    neighbor_offsets_[site_count_] = edge_count;
    neighbor_ids_.assign(edge_count, 0);
    neighbor_energies_.assign(edge_count, 0.0);

    double total_edge_energy = 0.0;
    for (int ix = 0; ix < config_.lattice.dimensions[0]; ++ix) {
        for (int iy = 0; iy < config_.lattice.dimensions[1]; ++iy) {
            for (int iz = 0; iz < config_.lattice.dimensions[2]; ++iz) {
                for (std::size_t basis = 0; basis < basis_count_; ++basis) {
                    const std::size_t site = site_id(ix, iy, iz, basis);
                    std::size_t edge = neighbor_offsets_[site];
                    for (const PatternNeighbor& neighbor : pattern[basis]) {
                        const int jx = wrap_index(ix + neighbor.dx, config_.lattice.dimensions[0]);
                        const int jy = wrap_index(iy + neighbor.dy, config_.lattice.dimensions[1]);
                        const int jz = wrap_index(iz + neighbor.dz, config_.lattice.dimensions[2]);
                        neighbor_ids_[edge] = static_cast<std::uint32_t>(
                            site_id(jx, jy, jz, neighbor.basis));
                        neighbor_energies_[edge] = neighbor.energy;
                        total_edge_energy += neighbor.energy;
                        ++edge;
                    }
                }
            }
        }
    }

    coordination_number_ = neighbor_offsets_[1] - neighbor_offsets_[0];
    cohesive_energy_ = 0.5 * total_edge_energy / static_cast<double>(site_count_);
}

void Simulation::initialize_state()
{
    std::fill(occupied_neighbor_count_.begin(), occupied_neighbor_count_.end(), 0);
    std::fill(local_field_.begin(), local_field_.end(), 0.0);

    for (std::size_t site = 0; site < site_count_; ++site) {
        if (occupied_[site]) {
            for (std::size_t edge = neighbor_offsets_[site]; edge < neighbor_offsets_[site + 1]; ++edge) {
                const std::size_t neighbor = neighbor_ids_[edge];
                ++occupied_neighbor_count_[neighbor];
                local_field_[neighbor] += neighbor_energies_[edge];
            }
        }
    }

    total_energy_ = 0.0;
    for (std::size_t site = 0; site < site_count_; ++site) {
        if (occupied_[site]) {
            total_energy_ += local_field_[site];
        }
    }
    total_energy_ *= 0.5;

    solid_interface_.reset(site_count_);
    solvent_interface_.reset(site_count_);
    evaporation_weights_.reset(site_count_);
    for (std::size_t site = 0; site < site_count_; ++site) {
        update_interface_membership(site);
    }
}

void Simulation::flip_site(std::size_t site_id)
{
    const bool was_occupied = occupied_[site_id] != 0;
    const int delta = was_occupied ? -1 : 1;
    total_energy_ += was_occupied ? -local_field_[site_id] : local_field_[site_id];
    occupied_[site_id] = static_cast<std::uint8_t>(!was_occupied);
    occupied_count_ = static_cast<std::size_t>(static_cast<long long>(occupied_count_) + delta);

    for (std::size_t edge = neighbor_offsets_[site_id]; edge < neighbor_offsets_[site_id + 1]; ++edge) {
        const std::size_t neighbor = neighbor_ids_[edge];
        occupied_neighbor_count_[neighbor] += delta;
        local_field_[neighbor] += static_cast<double>(delta) * neighbor_energies_[edge];
        update_interface_membership(neighbor);
    }
    update_interface_membership(site_id);
}

void Simulation::update_interface_membership(std::size_t site_id)
{
    const bool solid = occupied_[site_id]
        && (static_cast<int>(coordination_number_) - occupied_neighbor_count_[site_id] > 0);
    const bool solvent = !occupied_[site_id] && occupied_neighbor_count_[site_id] > 0;
    if (solid) {
        solid_interface_.add(site_id);
    } else {
        solid_interface_.remove(site_id);
    }
    if (solvent) {
        solvent_interface_.add(site_id);
    } else {
        solvent_interface_.remove(site_id);
    }
    evaporation_weights_.set(
        site_id,
        solid ? std::exp(config_.solvent.beta * local_field_[site_id]) : 0.0);
}

std::size_t Simulation::site_id(int ix, int iy, int iz, std::size_t basis) const
{
    const std::size_t ny = static_cast<std::size_t>(config_.lattice.dimensions[1]);
    const std::size_t nz = static_cast<std::size_t>(config_.lattice.dimensions[2]);
    return (((static_cast<std::size_t>(ix) * ny + static_cast<std::size_t>(iy)) * nz
             + static_cast<std::size_t>(iz))
            * basis_count_)
        + basis;
}

double Simulation::interaction_energy(double distance) const
{
    std::size_t interaction = 0;
    while (distance > config_.interactions.cutoffs[interaction]) {
        ++interaction;
    }
    return config_.interactions.interaction_energies[interaction];
}

} // namespace cgkmc
