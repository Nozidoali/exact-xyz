#include "prepare-state.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xyz {
namespace {

std::pair<uint32_t, uint32_t> maximize_difference_once(uint32_t num_qubits, std::unordered_set<uint32_t>& indices,
                                                       std::unordered_map<uint32_t, bool>& diff_values) {
    int                          max_diff = -1;
    std::unordered_set<uint32_t> max_diff_indices_1;
    uint32_t                     max_diff_qubit = 0;
    uint32_t                     max_diff_value = false;
    uint32_t                     length         = (uint32_t)indices.size();
    std::unordered_set<uint32_t> indices_1;

    for (uint32_t qubit = 0; qubit < num_qubits; qubit++) {
        if (diff_values.find(qubit) != diff_values.end())
            continue;
        indices_1.clear();
        for (auto index : indices)
            if ((index >> qubit) & 1u)
                indices_1.insert(index);

        int diff = abs((int)length - 2 * (int)indices_1.size());
        if (diff == (int)length)
            continue;
        if (diff > max_diff) {
            max_diff           = diff;
            max_diff_indices_1 = indices_1;
            max_diff_qubit     = qubit;
            max_diff_value     = length > 2 * indices_1.size();
        }
        if (max_diff == (int)length - 1)
            break;
    }

    if (max_diff_value)
        indices = max_diff_indices_1;
    else
        for (auto index : max_diff_indices_1)
            indices.erase(index);

    diff_values[max_diff_qubit] = max_diff_value;
    return {max_diff_qubit, max_diff_value};
}

} // namespace

ReductionResult cardinality_reduction_by_one(const QRState& state) {
    QRState                      new_state = state.clone();
    std::unordered_set<uint32_t> indices;
    for (auto [index, weight] : state.index_to_weight)
        indices.insert(index);

    std::unordered_map<uint32_t, bool> diff_values;
    uint32_t                           diff_qubit = 0, diff_value = 0;
    while (indices.size() > 1)
        std::tie(diff_qubit, diff_value) = maximize_difference_once(state.n_bits, indices, diff_values);

    uint32_t index0 = *indices.begin();
    diff_values.erase(diff_qubit);

    std::unordered_set<uint32_t> candidates;
    for (auto [index, weight] : state.index_to_weight) {
        if (indices.find(index) != indices.end())
            continue;
        bool valid = true;
        for (auto [qubit, value] : diff_values)
            if (((index >> qubit) & 1u) != (uint32_t)value) {
                valid = false;
                break;
            }
        if (valid)
            candidates.insert(index);
    }
    while (candidates.size() > 1)
        (void)maximize_difference_once(state.n_bits, candidates, diff_values);

    uint32_t                            index1 = *candidates.begin();
    std::vector<std::shared_ptr<QGate>> gates;

    for (uint32_t qubit = 0; qubit < state.n_bits; qubit++) {
        if (((index0 >> qubit) & 1u) == ((index1 >> qubit) & 1u))
            continue;
        if (qubit == diff_qubit)
            continue;
        auto gate = std::make_shared<CX>(diff_qubit, diff_value, qubit);
        gates.push_back(gate);
        new_state = (*gate)(new_state, true);
    }

    std::vector<uint32_t> ctrls;
    std::vector<bool>     phases;
    for (auto [qubit, value] : diff_values) {
        ctrls.push_back(qubit);
        phases.push_back(value);
    }

    uint32_t idx0 = index1 & (~(1u << diff_qubit));
    uint32_t idx1 = index1 | (1u << diff_qubit);
    auto     it0  = new_state.index_to_weight.find(idx0);
    auto     it1  = new_state.index_to_weight.find(idx1);
    assert(it0 != new_state.index_to_weight.end());
    assert(it1 != new_state.index_to_weight.end());

    double theta = 2.0 * atan2l((long double)it1->second, (long double)it0->second);
    if (index1 & (1u << diff_qubit))
        theta = -M_PI + theta;
    auto mcry_gate = std::make_shared<MCRY>(ctrls, phases, theta, diff_qubit);
    gates.push_back(mcry_gate);
    new_state = (*mcry_gate)(new_state, true);

    std::reverse(gates.begin(), gates.end());
    return {new_state, gates};
}

QCircuit prepare_sparse_state(const QRState& state) {
    QCircuit circuit(state.n_bits);
    QRState  curr_state = state.clone();
    while (curr_state.cardinality() > 1) {
        auto result = cardinality_reduction_by_one(curr_state);
        for (const auto& gate : result.gates)
            circuit.add_gate(gate);
        curr_state = result.state;
    }
    uint32_t index = curr_state.index_to_weight.begin()->first;
    for (uint32_t qubit = 0; qubit < state.n_bits; qubit++)
        if ((index >> qubit) & 1u)
            circuit.add_gate(std::make_shared<X>(qubit));
    circuit.reverse();
    return circuit;
}

} // namespace xyz
