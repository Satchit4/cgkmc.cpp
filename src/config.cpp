#include "cgkmc_core/config.hpp"

#include <cstdlib>
#include <fstream>
#include <istream>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>

namespace cgkmc {
namespace {

struct JsonValue;

using JsonArray = std::vector<JsonValue>;
using JsonObject = std::unordered_map<std::string, JsonValue>;

struct JsonValue {
    std::variant<double, bool, std::string, JsonArray, JsonObject> value;

    [[nodiscard]] const JsonArray& array() const
    {
        return std::get<JsonArray>(value);
    }

    [[nodiscard]] const JsonObject& object() const
    {
        return std::get<JsonObject>(value);
    }

    [[nodiscard]] double number() const
    {
        return std::get<double>(value);
    }
};

class JsonParser {
public:
    explicit JsonParser(std::string text)
        : text_(std::move(text))
    {
    }

    [[nodiscard]] JsonValue parse()
    {
        skip_ws();
        return parse_value();
    }

private:
    std::string text_;
    std::size_t pos_{};

    void skip_ws()
    {
        while (pos_ < text_.size() && static_cast<unsigned char>(text_[pos_]) <= ' ') {
            ++pos_;
        }
    }

    [[nodiscard]] JsonValue parse_value()
    {
        skip_ws();
        const char c = text_[pos_];
        if (c == '{') {
            return JsonValue{parse_object()};
        }
        if (c == '[') {
            return JsonValue{parse_array()};
        }
        if (c == '"') {
            return JsonValue{parse_string()};
        }
        if (c == 't') {
            pos_ += 4;
            return JsonValue{true};
        }
        if (c == 'f') {
            pos_ += 5;
            return JsonValue{false};
        }
        return JsonValue{parse_number()};
    }

    [[nodiscard]] JsonObject parse_object()
    {
        JsonObject object;
        ++pos_;
        skip_ws();
        while (text_[pos_] != '}') {
            const std::string key = parse_string();
            skip_ws();
            ++pos_;
            object.emplace(key, parse_value());
            skip_ws();
            if (text_[pos_] == ',') {
                ++pos_;
                skip_ws();
            }
        }
        ++pos_;
        return object;
    }

    [[nodiscard]] JsonArray parse_array()
    {
        JsonArray array;
        ++pos_;
        skip_ws();
        while (text_[pos_] != ']') {
            array.push_back(parse_value());
            skip_ws();
            if (text_[pos_] == ',') {
                ++pos_;
                skip_ws();
            }
        }
        ++pos_;
        return array;
    }

    [[nodiscard]] std::string parse_string()
    {
        std::string value;
        ++pos_;
        while (text_[pos_] != '"') {
            if (text_[pos_] == '\\') {
                ++pos_;
            }
            value.push_back(text_[pos_]);
            ++pos_;
        }
        ++pos_;
        return value;
    }

    [[nodiscard]] double parse_number()
    {
        const char* begin = text_.c_str() + pos_;
        char* end = nullptr;
        const double value = std::strtod(begin, &end);
        pos_ = static_cast<std::size_t>(end - text_.c_str());
        return value;
    }
};

[[nodiscard]] std::vector<double> read_number_vector(const JsonValue& value)
{
    std::vector<double> result;
    for (const JsonValue& element : value.array()) {
        result.push_back(element.number());
    }
    return result;
}

[[nodiscard]] std::array<int, 3> read_int3(const JsonValue& value)
{
    const JsonArray& array = value.array();
    return {
        static_cast<int>(array[0].number()),
        static_cast<int>(array[1].number()),
        static_cast<int>(array[2].number()),
    };
}

[[nodiscard]] std::array<double, 3> read_double3(const JsonValue& value)
{
    const JsonArray& array = value.array();
    return {array[0].number(), array[1].number(), array[2].number()};
}

[[nodiscard]] std::vector<std::array<double, 3>> read_basis(const JsonValue& value)
{
    std::vector<std::array<double, 3>> basis;
    for (const JsonValue& row : value.array()) {
        basis.push_back(read_double3(row));
    }
    return basis;
}

} // namespace

std::size_t LatticeConfig::basis_count() const
{
    return atomic_basis.size();
}

std::size_t LatticeConfig::site_count() const
{
    return static_cast<std::size_t>(dimensions[0])
        * static_cast<std::size_t>(dimensions[1])
        * static_cast<std::size_t>(dimensions[2])
        * basis_count();
}

double LatticeConfig::density() const
{
    return static_cast<double>(basis_count())
        / (lattice_parameters[0] * lattice_parameters[1] * lattice_parameters[2]);
}

double LatticeConfig::molecular_volume() const
{
    return 1.0 / density();
}

SimulationConfig load_config(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const JsonValue root = JsonParser(std::move(text)).parse();
    const JsonObject& object = root.object();
    const JsonObject& lattice = object.at("lattice").object();
    const JsonObject& interactions = object.at("interactions").object();
    const JsonObject& solvent = object.at("solvent").object();
    const JsonObject& growth = object.at("growth").object();

    SimulationConfig config;
    config.lattice.dimensions = read_int3(lattice.at("dimensions"));
    config.lattice.lattice_parameters = read_double3(lattice.at("lattice_parameters"));
    config.lattice.atomic_basis = read_basis(lattice.at("atomic_basis"));
    config.interactions.cutoffs = read_number_vector(interactions.at("cutoffs"));
    config.interactions.interaction_energies = read_number_vector(interactions.at("interaction_energies"));
    config.solvent.beta = solvent.at("beta").number();
    config.solvent.diffusivity = solvent.at("diffusivity").number();
    config.solvent.solubility_limit = solvent.at("solubility_limit").number();
    config.growth.initial_radius = growth.at("initial_radius").number();
    config.growth.num_steps = static_cast<std::size_t>(growth.at("num_steps").number());
    config.growth.desired_size = static_cast<std::size_t>(growth.at("desired_size").number());
    return config;
}

} // namespace cgkmc
