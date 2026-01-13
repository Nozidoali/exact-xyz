#include "qcircuit.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace xyz {
namespace {

std::vector<uint64_t> get_qubit_signatures(const QRState& state) {
    std::vector<uint64_t> signatures(state.n_bits, 0);
    for (const auto& [index, weight] : state.index_to_weight) {
        for (uint32_t j = 0; j < state.n_bits; j++) {
            signatures[j] = (signatures[j] << 1) | ((index >> j) & 1);
        }
    }
    return signatures;
}

uint64_t get_const1_signature(const QRState& state) {
    return (1ull << state.cardinality()) - 1;
}

std::optional<double> get_ap_ry_angles(const QRState& state, uint32_t qubit_index) {
    std::optional<double> theta;
    for (const auto& [index, weight] : state.index_to_weight) {
        uint32_t reversed_index = index ^ (1 << qubit_index);
        if (state.index_to_weight.find(reversed_index) == state.index_to_weight.end()) {
            return std::nullopt;
        }
        uint32_t index0 = index & ~(1 << qubit_index);
        uint32_t index1 = index | (1 << qubit_index);
        auto     it0    = state.index_to_weight.find(index0);
        auto     it1    = state.index_to_weight.find(index1);
        if (it0 == state.index_to_weight.end() || it1 == state.index_to_weight.end()) {
            return std::nullopt;
        }
        double weight0 = it0->second;
        double weight1 = it1->second;
        double _theta  = 2.0 * std::atan(weight1 / weight0);
        if (!theta.has_value()) {
            theta = _theta;
        } else if (std::abs(theta.value() - _theta) < 1e-10) {
            continue;
        } else {
            return std::nullopt;
        }
    }
    return theta;
}

QRState apply_ry_inverse(const QRState& state, uint32_t target, double theta) {
    double                     cos_half = std::cos(-theta / 2.0);
    double                     sin_half = std::sin(-theta / 2.0);
    std::map<uint32_t, double> temp_weights;
    for (const auto& [idx, _] : state.index_to_weight) {
        temp_weights[idx] = 0.0;
    }
    for (const auto& [idx, weight] : state.index_to_weight) {
        uint32_t rdx = idx ^ (1 << target);
        if (temp_weights.find(rdx) == temp_weights.end()) {
            temp_weights[rdx] = 0.0;
        }
        if (((idx >> target) & 1) == 0) {
            temp_weights[idx] += weight * cos_half;
            temp_weights[rdx] += weight * sin_half;
        } else {
            temp_weights[idx] += weight * cos_half;
            temp_weights[rdx] -= weight * sin_half;
        }
    }
    std::map<uint32_t, double> new_weights;
    for (const auto& [idx, w] : temp_weights) {
        if (std::abs(w) > QRState::eps) {
            new_weights[idx] = w;
        }
    }
    return QRState(new_weights, state.n_bits);
}

QRState apply_x(const QRState& state, uint32_t target) {
    std::map<uint32_t, double> new_weights;
    for (const auto& [index, weight] : state.index_to_weight) {
        new_weights[index ^ (1 << target)] = weight;
    }
    return QRState(new_weights, state.n_bits);
}

QRState apply_cx(const QRState& state, uint32_t control, bool control_phase, uint32_t target) {
    std::map<uint32_t, double> new_weights;
    for (const auto& [index, weight] : state.index_to_weight) {
        bool     ctrl_value = ((index >> control) & 1) != 0;
        uint32_t new_index  = (ctrl_value == control_phase) ? (index ^ (1 << target)) : index;
        new_weights[new_index] = weight;
    }
    return QRState(new_weights, state.n_bits);
}

struct SupportReductionResult {
    QRState                             state;
    std::vector<std::shared_ptr<QGate>> gates;
};

SupportReductionResult ry_reduction(const QRState& input_state) {
    QRState                             state = input_state;
    std::vector<std::shared_ptr<QGate>> gates;
    for (uint32_t qubit_index = 0; qubit_index < state.n_bits; qubit_index++) {
        auto theta = get_ap_ry_angles(state, qubit_index);
        if (theta.has_value()) {
            gates.push_back(std::make_shared<RY>(qubit_index, theta.value()));
            state = apply_ry_inverse(state, qubit_index, theta.value());
        }
    }
    return {state, gates};
}

SupportReductionResult x_reduction(const QRState& input_state, bool enable_cnot) {
    auto                                   signatures = get_qubit_signatures(input_state);
    auto                                   const1     = get_const1_signature(input_state);
    std::unordered_map<uint64_t, uint32_t> signature_to_qubits;
    QRState                                state = input_state;
    std::vector<std::shared_ptr<QGate>>    gates;
    for (uint32_t qubit_index = 0; qubit_index < signatures.size(); qubit_index++) {
        uint64_t signature = signatures[qubit_index];
        if (signature == 0) {
            continue;
        }
        if (signature == const1) {
            gates.push_back(std::make_shared<X>(qubit_index));
            state = apply_x(state, qubit_index);
            continue;
        }
        if (enable_cnot && signature_to_qubits.find(signature) != signature_to_qubits.end()) {
            uint32_t control_qubit = signature_to_qubits[signature];
            gates.push_back(std::make_shared<CX>(control_qubit, true, qubit_index));
            state = apply_cx(state, control_qubit, true, qubit_index);
            continue;
        }
        if (enable_cnot && signature_to_qubits.find(signature ^ const1) != signature_to_qubits.end()) {
            uint32_t control_qubit = signature_to_qubits[signature ^ const1];
            gates.push_back(std::make_shared<CX>(control_qubit, false, qubit_index));
            state = apply_cx(state, control_qubit, false, qubit_index);
            continue;
        }
        if (enable_cnot) {
            bool found = false;
            for (uint32_t q2 = qubit_index + 1; q2 < signatures.size(); q2++) {
                uint64_t sig2 = signatures[q2];
                if (signature_to_qubits.find(sig2 ^ signature) != signature_to_qubits.end()) {
                    uint32_t ctrl = signature_to_qubits[sig2 ^ signature];
                    gates.push_back(std::make_shared<CX>(ctrl, true, q2));
                    state = apply_cx(state, ctrl, true, q2);
                    gates.push_back(std::make_shared<CX>(q2, true, qubit_index));
                    state = apply_cx(state, q2, true, qubit_index);
                    found = true;
                    break;
                }
                if (signature_to_qubits.find(sig2 ^ const1 ^ signature) != signature_to_qubits.end()) {
                    uint32_t ctrl = signature_to_qubits[sig2 ^ const1 ^ signature];
                    gates.push_back(std::make_shared<CX>(ctrl, true, q2));
                    state = apply_cx(state, ctrl, true, q2);
                    gates.push_back(std::make_shared<CX>(q2, false, qubit_index));
                    state = apply_cx(state, q2, false, qubit_index);
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

SupportReductionResult support_reduction(const QRState& input_state) {
    auto x_result  = x_reduction(input_state, true);
    auto ry_result = ry_reduction(x_result.state);
    ry_result.gates.insert(ry_result.gates.end(), x_result.gates.begin(), x_result.gates.end());
    return ry_result;
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

} // namespace

QCircuit prepare_state_dense(const QRState& state) {
    auto     result = support_reduction(state);
    QCircuit circuit(state.n_bits);
    if (result.state.cardinality() == 1) {
        uint32_t index     = result.state.index_to_weight.begin()->first;
        bool     is_ground = (index == 0 && std::abs(result.state.index_to_weight.at(0) - 1.0) < 1e-6);
        if (is_ground) {
            for (const auto& gate : result.gates) {
                circuit.add_gate(gate);
            }
        } else {
            for (uint32_t target = 0; target < result.state.n_bits; target++) {
                if ((index >> target) & 1) {
                    circuit.add_gate(std::make_shared<X>(target));
                }
            }
            for (auto it = result.gates.rbegin(); it != result.gates.rend(); ++it) {
                circuit.add_gate(conjugate_gate(*it));
            }
        }
    } else if (result.state.cardinality() > 1) {
        auto sparse_circuit = prepare_sparse_state(result.state);
        for (const auto& gate : sparse_circuit.pGates) {
            circuit.add_gate(gate);
        }
        for (auto it = result.gates.rbegin(); it != result.gates.rend(); ++it) {
            circuit.add_gate(conjugate_gate(*it));
        }
    }
    return circuit;
}

} // namespace xyz
