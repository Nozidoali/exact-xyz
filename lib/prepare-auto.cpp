#include "prepare-state.hpp"

#include <iostream>
#include <memory>
#include <vector>

namespace xyz {
namespace {

bool is_uniform_state(const QRState& state) {
    if (state.index_to_weight.empty())
        return true;
    double reference_weight = state.index_to_weight.begin()->second;
    for (const auto& [index, weight] : state.index_to_weight) {
        if (std::abs(weight - reference_weight) > QRState::eps)
            return false;
    }
    return true;
}

void prepare_auto_rec(const QRState& state, std::vector<std::shared_ptr<QGate>>& gates, bool verbose) {
    std::vector<std::shared_ptr<QGate>> support_reduction_gates;
    QRState                             reduced_state = state;

    auto support_result     = support_reduction(state);
    reduced_state           = support_result.state;
    support_reduction_gates = support_result.gates;

    std::vector<uint32_t> supports     = reduced_state.get_supports();
    uint32_t              num_supports = supports.size();
    uint32_t              cardinality  = reduced_state.cardinality();

    if (verbose)
        std::cout << "n=" << num_supports << " card=" << cardinality << "\n";

    if (cardinality == 1) {
        uint32_t index = reduced_state.index_to_weight.begin()->first;
        for (uint32_t qubit = 0; qubit < reduced_state.n_bits; qubit++)
            if ((index >> qubit) & 1)
                support_reduction_gates.push_back(std::make_shared<X>(qubit));
        gates.insert(gates.end(), support_reduction_gates.begin(), support_reduction_gates.end());
        return;
    }

    const uint32_t EXACT_SYNTHESIS_DENSITY_THRESHOLD = 100;
    const uint32_t EXACT_SYNTHESIS_CNOT_LIMIT        = 100;
    const uint32_t n_qubits_max                      = 5;

    if (num_supports <= n_qubits_max && cardinality <= EXACT_SYNTHESIS_DENSITY_THRESHOLD) {
        QCircuit   circ(reduced_state.n_bits);
        bfs_params params;
        bool       bfs_success = prepare_state_bfs(reduced_state, circ, params, false);

        if (bfs_success) {
            std::vector<std::shared_ptr<QGate>> exact_gates;
            for (auto it = circ.pGates.rbegin(); it != circ.pGates.rend(); ++it)
                exact_gates.push_back(*it);
            gates.insert(gates.end(), exact_gates.begin(), exact_gates.end());
            gates.insert(gates.end(), support_reduction_gates.begin(), support_reduction_gates.end());
            return;
        }
    }

    std::vector<std::shared_ptr<QGate>> cardinality_reduction_gates;
    auto                                card_result = cardinality_reduction_by_one(reduced_state);
    cardinality_reduction_gates                     = card_result.gates;
    QRState card_reduced_state                      = card_result.state;

    std::vector<std::shared_ptr<QGate>> rec_gates;
    prepare_auto_rec(card_reduced_state, rec_gates, verbose);

    gates.insert(gates.end(), rec_gates.begin(), rec_gates.end());
    gates.insert(gates.end(), cardinality_reduction_gates.begin(), cardinality_reduction_gates.end());
    gates.insert(gates.end(), support_reduction_gates.begin(), support_reduction_gates.end());
}

} // namespace

QCircuit prepare_state_auto(const QRState& state, bool verbose) {
    QCircuit                            circuit(state.n_bits);
    std::vector<std::shared_ptr<QGate>> gates;
    prepare_auto_rec(state, gates, verbose);
    for (auto it = gates.rbegin(); it != gates.rend(); ++it)
        circuit.add_gate(*it);
    return circuit;
}

} // namespace xyz
