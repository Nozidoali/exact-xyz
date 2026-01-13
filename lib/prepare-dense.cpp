#include "prepare-state.hpp"

#include <cassert>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xyz {
namespace {

ReductionResult ry_reduction(const QRState& input_state) {
    QRState                             state = input_state;
    std::vector<std::shared_ptr<QGate>> gates;
    for (uint32_t qubit_index = 0; qubit_index < state.n_bits; qubit_index++) {
        auto theta = state.get_ap_ry_angles(qubit_index);
        if (theta.has_value()) {
            auto ry_gate = std::make_shared<RY>(qubit_index, theta.value());
            gates.push_back(ry_gate);
            state = (*ry_gate)(state, true);
        }
    }
    return {state, gates};
}

ReductionResult x_reduction(const QRState& input_state, bool enable_cnot) {
    auto                                   signatures = input_state.get_qubit_signatures();
    auto                                   const1     = input_state.get_const1_signature();
    std::unordered_map<uint64_t, uint32_t> signature_to_qubits;
    QRState                                state = input_state;
    std::vector<std::shared_ptr<QGate>>    gates;
    for (uint32_t qubit_index = 0; qubit_index < signatures.size(); qubit_index++) {
        uint64_t signature = signatures[qubit_index];
        if (signature == 0) {
            continue;
        }
        if (signature == const1) {
            auto x_gate = std::make_shared<X>(qubit_index);
            gates.push_back(x_gate);
            state = (*x_gate)(state, false);
            continue;
        }
        if (enable_cnot && signature_to_qubits.find(signature) != signature_to_qubits.end()) {
            uint32_t control_qubit = signature_to_qubits[signature];
            auto     cx_gate       = std::make_shared<CX>(control_qubit, true, qubit_index);
            gates.push_back(cx_gate);
            state = (*cx_gate)(state, false);
            continue;
        }
        if (enable_cnot && signature_to_qubits.find(signature ^ const1) != signature_to_qubits.end()) {
            uint32_t control_qubit = signature_to_qubits[signature ^ const1];
            auto     cx_gate       = std::make_shared<CX>(control_qubit, false, qubit_index);
            gates.push_back(cx_gate);
            state = (*cx_gate)(state, false);
            continue;
        }
        if (enable_cnot) {
            bool found = false;
            for (uint32_t q2 = qubit_index + 1; q2 < signatures.size(); q2++) {
                uint64_t sig2 = signatures[q2];
                if (signature_to_qubits.find(sig2 ^ signature) != signature_to_qubits.end()) {
                    uint32_t ctrl = signature_to_qubits[sig2 ^ signature];
                    auto     cx1  = std::make_shared<CX>(ctrl, true, q2);
                    gates.push_back(cx1);
                    state    = (*cx1)(state, false);
                    auto cx2 = std::make_shared<CX>(q2, true, qubit_index);
                    gates.push_back(cx2);
                    state = (*cx2)(state, false);
                    found = true;
                    break;
                }
                if (signature_to_qubits.find(sig2 ^ const1 ^ signature) != signature_to_qubits.end()) {
                    uint32_t ctrl = signature_to_qubits[sig2 ^ const1 ^ signature];
                    auto     cx1  = std::make_shared<CX>(ctrl, true, q2);
                    gates.push_back(cx1);
                    state    = (*cx1)(state, false);
                    auto cx2 = std::make_shared<CX>(q2, false, qubit_index);
                    gates.push_back(cx2);
                    state = (*cx2)(state, false);
                    found = true;
                    break;
                }
            }
            if (found) {
                continue;
            }
        }
        signature_to_qubits[signature] = qubit_index;
    }
    return {state, gates};
}

std::shared_ptr<QGate> conjugate_gate(const std::shared_ptr<QGate>& gate) {
    if (auto ry = std::dynamic_pointer_cast<RY>(gate)) {
        return std::make_shared<RY>(ry->target, -ry->theta);
    } else if (auto cry = std::dynamic_pointer_cast<CRY>(gate)) {
        return std::make_shared<CRY>(cry->ctrl, cry->phase, -cry->theta, cry->target);
    } else if (auto mcry = std::dynamic_pointer_cast<MCRY>(gate)) {
        return std::make_shared<MCRY>(mcry->ctrls, mcry->phases, -mcry->theta, mcry->target);
    }
    return gate;
}

uint32_t select_informative_qubit(const QRState& state, const std::vector<uint32_t>& supports) {
    uint32_t best_qubit     = 0;
    uint32_t min_difference = UINT32_MAX;
    uint32_t length         = state.cardinality();

    assert(length >= 2);

    for (uint32_t qubit : supports) {
        uint32_t length0 = 0;
        for (const auto& [index, weight] : state.index_to_weight) {
            if ((index & (1u << qubit)) == 0) {
                length0++;
            }
        }

        uint32_t difference = abs((int)length - 2 * (int)length0);
        if (difference < min_difference) {
            min_difference = difference;
            best_qubit     = qubit;
        }
    }

    assert(best_qubit < state.n_bits);
    return best_qubit;
}

std::pair<std::vector<double>, std::vector<uint32_t>>
rotation_angles_optimization(const std::vector<double>& rotation_angles, const std::vector<uint32_t>& control_indices) {

    assert(rotation_angles.size() == (1u << control_indices.size()));

    constexpr double             MIN_ROTATION_ANGLE_SEPARATION = 1e-6;
    std::unordered_set<uint32_t> dont_cares;

    for (uint32_t index = 0; index < control_indices.size(); index++) {
        bool is_dont_care = true;
        for (uint32_t rotation_index = 0; rotation_index < rotation_angles.size(); rotation_index++) {
            uint32_t reversed_index = rotation_index ^ (1u << index);
            if (std::abs(rotation_angles[reversed_index] - rotation_angles[rotation_index]) >
                MIN_ROTATION_ANGLE_SEPARATION) {
                is_dont_care = false;
                break;
            }
        }
        if (is_dont_care) {
            dont_cares.insert(index);
        }
    }

    if (dont_cares.empty()) {
        return {rotation_angles, control_indices};
    }

    std::vector<uint32_t> new_control_indices;
    std::vector<uint32_t> old_indices;
    for (uint32_t old_index = 0; old_index < control_indices.size(); old_index++) {
        if (dont_cares.find(old_index) == dont_cares.end()) {
            new_control_indices.push_back(control_indices[old_index]);
            old_indices.push_back(old_index);
        }
    }

    std::vector<double> new_rotation_angles;
    for (uint32_t new_index = 0; new_index < (1u << new_control_indices.size()); new_index++) {
        uint32_t old_index = 0;
        for (uint32_t i = 0; i < new_control_indices.size(); i++) {
            old_index |= ((new_index >> i) & 1u) << old_indices[i];
        }
        new_rotation_angles.push_back(rotation_angles[old_index]);
    }

    return {new_rotation_angles, new_control_indices};
}

ReductionResult qubit_reduction_by_one(const QRState& state) {
    auto supports = state.get_supports();
    if (supports.size() <= 1) {
        return {state, {}};
    }

    uint32_t pivot = select_informative_qubit(state, supports);

    std::vector<uint32_t> control_indices;
    for (uint32_t support : supports) {
        if (support != pivot) {
            control_indices.push_back(support);
        }
    }

    std::vector<std::pair<double, double>> rotation_table(1u << control_indices.size(), {0.0, 0.0});

    for (const auto& [index, weight] : state.index_to_weight) {
        uint32_t rotation_index = 0;
        for (uint32_t i = 0; i < control_indices.size(); i++) {
            if (index & (1u << control_indices[i])) {
                rotation_index |= 1u << i;
            }
        }

        if (index & (1u << pivot)) {
            rotation_table[rotation_index].second += weight;
        } else {
            rotation_table[rotation_index].first += weight;
        }
    }

    std::vector<double> rotation_angles;
    for (const auto& entry : rotation_table) {
        if (entry.first == 0.0) {
            rotation_angles.push_back(M_PI);
        } else if (entry.second == 0.0) {
            rotation_angles.push_back(0.0);
        } else {
            rotation_angles.push_back(2.0 * std::atan2(entry.second, entry.first));
        }
    }

    auto [optimized_angles, optimized_controls] = rotation_angles_optimization(rotation_angles, control_indices);

    std::vector<std::shared_ptr<QGate>> gates;
    std::shared_ptr<QGate>              gate;

    if (optimized_angles.size() == 1) {
        gate = std::make_shared<RY>(pivot, optimized_angles[0]);
    } else {
        std::vector<bool> phases(optimized_controls.size(), true);
        gate = std::make_shared<MCMY>(optimized_controls, phases, optimized_angles, pivot);
    }
    gates.push_back(gate);

    std::map<uint32_t, double>   new_weights;
    std::unordered_set<uint32_t> processed;

    for (const auto& [index, weight] : state.index_to_weight) {
        if (processed.find(index) != processed.end()) {
            continue;
        }

        uint32_t idx0         = index & ~(1u << pivot);
        uint32_t idx1         = index | (1u << pivot);
        uint32_t reversed_idx = index ^ (1u << pivot);

        if (state.index_to_weight.find(reversed_idx) == state.index_to_weight.end()) {
            new_weights[idx0] = weight;
        } else {
            double weight0       = state.index_to_weight.at(idx0);
            double weight1       = state.index_to_weight.at(idx1);
            double merged_weight = std::sqrt(weight0 * weight0 + weight1 * weight1);
            new_weights[idx0]    = merged_weight;
            processed.insert(idx1);
        }
        processed.insert(index);
    }

    for (auto it = new_weights.begin(); it != new_weights.end();) {
        if (std::abs(it->second) < QRState::eps) {
            it = new_weights.erase(it);
        } else {
            ++it;
        }
    }

    QRState reduced_state(new_weights, state.n_bits);

    return {reduced_state, gates};
}

} // namespace

QCircuit prepare_state_dense(const QRState& state) {
    QCircuit circuit(state.n_bits);
    QRState  curr_state = state.clone();

    std::vector<uint32_t> supports = curr_state.get_supports();
    while (supports.size() > 1) {
        auto result = qubit_reduction_by_one(curr_state);
        for (const auto& gate : result.gates) {
            circuit.add_gate(gate);
        }
        curr_state = result.state;
        supports   = curr_state.get_supports();
    }

    if (curr_state.cardinality() == 1) {
        uint32_t index     = curr_state.index_to_weight.begin()->first;
        bool     is_ground = (index == 0 && std::abs(curr_state.index_to_weight.at(0) - 1.0) < 1e-6);
        if (!is_ground) {
            for (uint32_t target = 0; target < curr_state.n_bits; target++) {
                if ((index >> target) & 1) {
                    circuit.add_gate(std::make_shared<X>(target));
                }
            }
        }
    } else if (curr_state.cardinality() > 1) {
        auto sparse_circuit = prepare_sparse_state(curr_state);
        for (const auto& gate : sparse_circuit.pGates) {
            circuit.add_gate(gate);
        }
    }

    circuit.reverse();
    return circuit;
}

ReductionResult support_reduction(const QRState& input_state) {
    auto x_result  = x_reduction(input_state, true);
    auto ry_result = ry_reduction(x_result.state);
    ry_result.gates.insert(ry_result.gates.end(), x_result.gates.begin(), x_result.gates.end());
    return ry_result;
}

} // namespace xyz
