#include "qcircuit.hpp"

#include <cassert>
#include <cmath>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xyz {
namespace {
template <bool verbose>
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

template <bool verbose> QRState cardinality_reduction_impl(QCircuit& qcircuit, const QRState& state) {
    QRState                      new_state = state.clone();
    std::unordered_set<uint32_t> indices;
    for (auto [index, weight] : state.index_to_weight)
        indices.insert(index);

    std::unordered_map<uint32_t, bool> diff_values;
    uint32_t                           diff_qubit = 0, diff_value = 0;
    while (indices.size() > 1)
        std::tie(diff_qubit, diff_value) = maximize_difference_once<verbose>(state.n_bits, indices, diff_values);

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
        (void)maximize_difference_once<verbose>(state.n_bits, candidates, diff_values);

    uint32_t index1 = *candidates.begin();
    for (uint32_t qubit = 0; qubit < state.n_bits; qubit++) {
        if (((index0 >> qubit) & 1u) == ((index1 >> qubit) & 1u))
            continue;
        if (qubit == diff_qubit)
            continue;
        qcircuit.add_gate(std::make_shared<CX>(diff_qubit, diff_value, qubit));
        new_state = (*std::make_shared<CX>(diff_qubit, diff_value, qubit))(new_state, true);
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
    qcircuit.add_gate(std::make_shared<MCRY>(ctrls, phases, theta, diff_qubit));
    new_state = (*std::make_shared<MCRY>(ctrls, phases, theta, diff_qubit))(new_state, true);
    return new_state;
}

void add_cx(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t control, uint32_t target) {
    if (qVal[control].has_value() && !qVal[control].value())
        return;
    if (qVal[control].has_value() && qVal[control].value()) {
        circuit.add_gate(std::make_shared<X>(target));
        if (qVal[target].has_value())
            qVal[target] = !qVal[target].value();
    } else {
        circuit.add_gate(std::make_shared<CX>(control, true, target));
        qVal[target] = std::nullopt;
    }
}

void add_ry(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t target, double theta) {
    circuit.add_gate(std::make_shared<RY>(target, theta));
    qVal[target] = std::nullopt;
}

void add_cry(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t control, uint32_t target,
             double theta) {
    if (qVal[control].has_value() && !qVal[control].value())
        return;
    if (qVal[control].has_value() && qVal[control].value())
        add_ry(circuit, qVal, target, theta);
    else {
        circuit.add_gate(std::make_shared<CRY>(control, true, theta, target));
        qVal[target] = std::nullopt;
    }
}

void add_mcry(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t control1, uint32_t control2,
              uint32_t target, double theta) {
    if (qVal[control1].has_value() && !qVal[control1].value())
        return;
    if (qVal[control2].has_value() && !qVal[control2].value())
        return;
    if (qVal[control1].has_value() && qVal[control1].value())
        return add_cry(circuit, qVal, control2, target, theta);
    if (qVal[control2].has_value() && qVal[control2].value())
        return add_cry(circuit, qVal, control1, target, theta);
    std::vector<uint32_t> controls = {control1, control2};
    circuit.add_gate(std::make_shared<MCRY>(controls, theta, target));
    qVal[target] = std::nullopt;
}

void insert_mu(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t n, uint32_t j) {
    double theta = 2 * std::acos(std::sqrt(1.0 / (n - j)));
    add_cx(circuit, qVal, j + 1, j);
    add_cry(circuit, qVal, j, j + 1, theta);
    add_cx(circuit, qVal, j + 1, j);
}

void insert_M(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t n, uint32_t k, uint32_t j,
              uint32_t i) {
    double theta = 2 * std::acos(std::sqrt((i + 1.0) / (n - j)));
    add_cx(circuit, qVal, j + i + 1, j);
    add_mcry(circuit, qVal, j + i, j, j + i + 1, theta);
    add_cx(circuit, qVal, j + i + 1, j);
}

void insert_scs(QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t n, uint32_t k, uint32_t j) {
    assert(k >= 1);
    insert_mu(circuit, qVal, n, j);
    for (uint32_t i = 1; i < k; ++i) {
        if (j + i + 1 >= n)
            break;
        insert_M(circuit, qVal, n, k, j, i);
    }
}
} // namespace

QCircuit prepare_sparse_state(const QRState& state) {
    QCircuit circuit(state.n_bits);
    QRState  curr_state = state.clone();
    while (curr_state.cardinality() > 1)
        curr_state = cardinality_reduction_impl<false>(circuit, curr_state);
    uint32_t index = curr_state.index_to_weight.begin()->first;
    for (uint32_t qubit = 0; qubit < state.n_bits; qubit++)
        if ((index >> qubit) & 1u)
            circuit.add_gate(std::make_shared<X>(qubit));
    circuit.reverse();
    return circuit;
}

QCircuit prepare_dicke_state(int n, int k) {
    QCircuit                         circuit(n);
    std::vector<std::optional<bool>> qVal(n, false);
    for (int i = 0; i < k; ++i) {
        circuit.add_gate(std::make_shared<X>(i));
        qVal[i] = true;
    }
    for (int i = 0; i < n - 1; ++i)
        insert_scs(circuit, qVal, (uint32_t)n, (uint32_t)k, (uint32_t)i);
    return circuit;
}

} // namespace xyz
