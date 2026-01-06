#include "qgate.hpp"

#include "qstate.hpp"

#include <cmath>
#include <stdexcept>

namespace xyz {
bool Rotation::is_trivial(double theta, bool use_x) {
    bool is_zero = std::abs(theta) < eps || std::abs(theta - 2 * M_PI) < eps;
    bool is_pi   = std::abs(theta - M_PI) < eps || std::abs(theta + M_PI) < eps;
    if (use_x)
        return is_zero || is_pi;
    return is_zero;
}

QRState X::operator()(const QRState& state, const bool reverse) const {
    (void)reverse;
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        uint32_t new_index                   = index ^ (1 << target);
        new_state.index_to_weight[new_index] = weight;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState RU2::operator()(const QRState& state, uint32_t target, const bool reverse) const {
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        if (new_state.index_to_weight.find(index) == new_state.index_to_weight.end())
            new_state.index_to_weight[index] = 0;
        uint32_t new_index = index ^ (1 << target);
        if (index & (1 << target)) {
            new_state.index_to_weight[new_index] += c10[reverse] * weight;
            new_state.index_to_weight[index] += c11[reverse] * weight;
        } else {
            new_state.index_to_weight[index] += c00[reverse] * weight;
            new_state.index_to_weight[new_index] += c01[reverse] * weight;
        }
    }
    for (auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end();) {
        if (std::abs(it->second) < QRState::eps)
            it = new_state.index_to_weight.erase(it);
        else
            ++it;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QState U2::operator()(const QState& state, uint32_t target, const bool reverse) const {
    QState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        if (new_state.index_to_weight.find(index) == new_state.index_to_weight.end())
            new_state.index_to_weight[index] = 0;
        uint32_t new_index = index ^ (1 << target);
        if (index & (1 << target)) {
            new_state.index_to_weight[new_index] += c10[reverse] * weight;
            new_state.index_to_weight[index] += c11[reverse] * weight;
        } else {
            new_state.index_to_weight[index] += c00[reverse] * weight;
            new_state.index_to_weight[new_index] += c01[reverse] * weight;
        }
    }
    for (auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end();) {
        if (std::abs(it->second) < QState::eps)
            it = new_state.index_to_weight.erase(it);
        else
            ++it;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState U2::operator()(const QRState& state, uint32_t target, const bool reverse) const {
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        if (new_state.index_to_weight.find(index) == new_state.index_to_weight.end())
            new_state.index_to_weight[index] = 0;
        uint32_t new_index = index ^ (1 << target);
        if (index & (1 << target)) {
            new_state.index_to_weight[new_index] += c10[reverse].real() * weight;
            new_state.index_to_weight[index] += c11[reverse].real() * weight;
        } else {
            new_state.index_to_weight[index] += c00[reverse].real() * weight;
            new_state.index_to_weight[new_index] += c01[reverse].real() * weight;
        }
    }
    for (auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end();) {
        if (std::abs(it->second) < QRState::eps)
            it = new_state.index_to_weight.erase(it);
        else
            ++it;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState CX::operator()(const QRState& state, const bool reverse) const {
    (void)reverse; // the conjugate of CX is CX
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        uint32_t new_index = index;
        if ((bool)((index >> ctrl) & (uint32_t)1u) == phase)
            new_index ^= (1 << target);
        new_state.index_to_weight[new_index] = weight;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState CCX::operator()(const QRState& state, const bool reverse) const {
    (void)reverse; // the conjugate of CCX is CCX
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        uint32_t new_index = index;
        if ((bool)((index >> ctrls[0]) & (uint32_t)1u) == phases[0] &&
            (bool)((index >> ctrls[1]) & (uint32_t)1u) == phases[1])
            new_index ^= (1 << target);
        new_state.index_to_weight[new_index] = weight;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState CRY::operator()(const QRState& state, const bool reverse) const {
    auto    _theta = reverse ? -theta : theta;
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        if ((bool)((index >> ctrl) & (uint32_t)1u) != phase) {
            new_state.index_to_weight[index] = weight;
            continue;
        }
        if (new_state.index_to_weight.find(index) == new_state.index_to_weight.end())
            new_state.index_to_weight[index] = 0;
        uint32_t new_index = index ^ (1 << target);
        if (index & (1 << target)) {
            new_state.index_to_weight[new_index] += c10[reverse] * weight;
            new_state.index_to_weight[index] += c11[reverse] * weight;
        } else {
            new_state.index_to_weight[index] += c00[reverse] * weight;
            new_state.index_to_weight[new_index] += c01[reverse] * weight;
        }
    }
    for (auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end();) {
        if (std::abs(it->second) < QRState::eps)
            it = new_state.index_to_weight.erase(it);
        else
            ++it;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState MCRY::operator()(const QRState& state, const bool reverse) const {
    QRState new_state;
    for (const auto& [index, weight] : state.index_to_weight) {
        bool valid = true;
        for (uint32_t i = 0; i < ctrls.size(); i++) {
            if ((bool)((index >> ctrls[i]) & (uint32_t)1) != phases[i]) {
                valid = false;
                break;
            }
        }
        if (!valid) {
            new_state.index_to_weight[index] = weight;
            continue;
        }
        if (new_state.index_to_weight.find(index) == new_state.index_to_weight.end())
            new_state.index_to_weight[index] = 0;
        uint32_t new_index = index ^ (1 << target);
        if (index & (1 << target)) {
            new_state.index_to_weight[new_index] += c10[reverse] * weight;
            new_state.index_to_weight[index] += c11[reverse] * weight;
        } else {
            new_state.index_to_weight[index] += c00[reverse] * weight;
            new_state.index_to_weight[new_index] += c01[reverse] * weight;
        }
    }
    for (auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end();) {
        if (std::abs(it->second) < QRState::eps)
            it = new_state.index_to_weight.erase(it);
        else
            ++it;
    }
    new_state.n_bits = state.n_bits;
    return new_state;
}

QRState S::operator()(const QRState& state, const bool reverse) const {
    (void)reverse;
    throw std::runtime_error("S gate not supported in real-amplitude QRState simulation");
}

QRState QROM_MCRY::operator()(const QRState& state, const bool reverse) const {
    QRState result = state.clone();
    for (uint32_t i = 0; i < rotation_table.size(); i++) {
        std::vector<uint32_t> active_ctrls;
        std::vector<bool>     active_phases;
        for (uint32_t j = 0; j < ctrls.size(); j++) {
            active_ctrls.push_back(ctrls[j]);
            active_phases.push_back((i >> j) & 1);
        }
        MCRY mcry(active_ctrls, active_phases, rotation_table[i], target);
        result = mcry(result, reverse);
    }
    return result;
}

} // namespace xyz