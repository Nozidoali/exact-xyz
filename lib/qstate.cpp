#include "qstate.hpp"

#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <unordered_set>

namespace xyz {
bool QRState::is_ground() const {
    return index_to_weight.size() == 1 && index_to_weight.find(0) != index_to_weight.end();
}
std::string QRState::to_string() const {
    std::string str;
    for (auto it = index_to_weight.begin(); it != index_to_weight.end(); it++) {
        auto index  = it->first;
        auto weight = it->second;
        str += std::to_string(weight) + "*|";
        for (uint32_t i = 0; i < n_bits; i++)
            str += ((index >> i) & 1) ? "1" : "0";
        str += ">";
        if (std::next(it) != index_to_weight.end())
            str += " + ";
    }
    return str;
}
std::ostream& operator<<(std::ostream& os, const QRState& obj) {
    os << obj.to_string();
    return os;
}
std::unordered_map<uint32_t, double> QRState::to_ry_table(uint32_t target) const {
    std::unordered_map<uint32_t, double> ry_table;
    for (const auto& [index, weight] : index_to_weight) {
        auto index_0 = index & (~(1 << target));
        auto index_1 = index | (1 << target);

        if (ry_table.find(index_0) != ry_table.end())
            continue;
        if (index_to_weight.find(index_0) == index_to_weight.end())
            ry_table[index_0] = M_PI;
        else if (index_to_weight.find(index_1) == index_to_weight.end())
            ry_table[index_0] = 0;
        else
            ry_table[index_0] = 2 * atan2(index_to_weight.at(index_1), index_to_weight.at(index_0));
    }
    return ry_table;
}
uint64_t QRState::repr() const {
    if (hash_value.has_value())
        return hash_value.value();
    uint64_t h   = 1469598103934665603ull;
    auto     mix = [&](uint64_t x) {
        h ^= x;
        h *= 1099511628211ull;
    };
    mix(n_bits);
    for (const auto& [index, weight] : index_to_weight) {
        mix(index);
        int64_t q = (int64_t)std::llround(weight / QRState::eps);
        mix((uint64_t)q);
    }
    hash_value = h;
    return h;
}
QRState QRState::clone() const {
    std::map<uint32_t, double> index_to_weight_copy;
    for (const auto& [index, weight] : index_to_weight)
        index_to_weight_copy[index] = weight;
    return QRState(index_to_weight_copy, n_bits);
}
QRState ground_rstate(uint32_t n_bits) {
    std::map<uint32_t, double> index_to_weight;
    index_to_weight[0] = 1.0;
    return QRState(index_to_weight, n_bits);
}
QRState dicke_state(uint32_t n, uint32_t k) {
    std::map<uint32_t, double> index_to_weight;
    double                     total_weight = 0;
    for (uint32_t i = 0; i < (1 << n); i++) {
        uint32_t count = 0;
        for (uint32_t j = 0; j < n; j++)
            count += (i >> j) & 1;
        if (count == k) {
            index_to_weight[i] = 1.0;
            total_weight += 1.0;
        }
    }
    for (auto& [index, weight] : index_to_weight)
        weight /= std::sqrt(total_weight);
    return QRState(index_to_weight, n);
}

QRState random_rstate(uint32_t n_bits, uint32_t cardinality, uint64_t seed) {
    if (n_bits >= 32)
        throw std::invalid_argument("random_rstate: n_bits must be < 32 (index type is uint32_t)");
    if (cardinality == 0)
        throw std::invalid_argument("random_rstate: cardinality must be >= 1");
    const uint64_t dim = 1ull << n_bits;
    if ((uint64_t)cardinality > dim)
        throw std::invalid_argument("random_rstate: cardinality exceeds Hilbert space dimension");

    uint64_t actual_seed = seed;
    if (actual_seed == 0) {
        std::random_device rd;
        actual_seed = ((uint64_t)rd() << 32) ^ (uint64_t)rd();
    }
    std::mt19937_64                         rng(actual_seed ? actual_seed : 1);
    std::uniform_int_distribution<uint32_t> pick_index(0u, (uint32_t)(dim - 1));
    std::normal_distribution<double>        pick_weight(0.0, 1.0);

    std::unordered_set<uint32_t> chosen;
    while (chosen.size() < cardinality) {
        chosen.insert(pick_index(rng));
    }

    std::map<uint32_t, double> index_to_weight;
    double                     norm2 = 0.0;
    for (uint32_t index : chosen) {
        double w               = pick_weight(rng);
        index_to_weight[index] = w;
        norm2 += w * w;
    }

    if (norm2 <= 0.0) {
        auto it = index_to_weight.begin();
        for (auto jt = std::next(it); jt != index_to_weight.end(); ++jt)
            jt->second = 0.0;
        if (it != index_to_weight.end())
            it->second = 1.0;
    } else {
        const double inv_norm = 1.0 / std::sqrt(norm2);
        for (auto& [index, w] : index_to_weight)
            w *= inv_norm;
    }

    return QRState(index_to_weight, n_bits);
}

std::ostream& operator<<(std::ostream& os, const QState& obj) {
    for (auto it = obj.index_to_weight.begin(); it != obj.index_to_weight.end(); it++) {
        auto index  = it->first;
        auto weight = it->second;
        os << weight << "*|";
        for (uint32_t i = 0; i < obj.n_bits; i++)
            os << ((index >> i) & 1);
        os << ">";
        if (std::next(it) != obj.index_to_weight.end())
            os << " + ";
    }
    return os;
}

QState ground_state(uint32_t n_bits) {
    std::map<uint32_t, std::complex<double>> index_to_weight;
    index_to_weight[0] = 1.0;
    return QState(index_to_weight, n_bits);
}
} // namespace xyz