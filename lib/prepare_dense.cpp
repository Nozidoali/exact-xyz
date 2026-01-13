#include "qcircuit.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace xyz {
namespace {

uint32_t select_pivot_qubit(const QRState& state, const std::vector<uint32_t>& supports) {
    int      max_diff   = -1;
    uint32_t best_qubit = supports[0];
    uint32_t length     = state.cardinality();

    for (uint32_t qubit : supports) {
        uint32_t count_0 = 0;
        for (const auto& [index, weight] : state.index_to_weight)
            if ((index & (1 << qubit)) == 0)
                count_0++;

        int diff = std::abs((int)length - 2 * (int)count_0);
        if (diff > max_diff) {
            max_diff   = diff;
            best_qubit = qubit;
        }
        if (diff == (int)length - 1)
            break;
    }

    return best_qubit;
}

struct Cofactors {
    QRState neg_state;
    QRState pos_state;
    double  weight0;
    double  weight1;
};

Cofactors compute_cofactors(const QRState& state, uint32_t pivot) {
    std::map<uint32_t, double> neg_map, pos_map;
    double                     w0 = 0.0, w1 = 0.0;

    for (const auto& [index, weight] : state.index_to_weight) {
        uint32_t idx_0 = index & ~(1 << pivot);
        if (index & (1 << pivot)) {
            pos_map[idx_0] = weight;
            w1 += weight * weight;
        } else {
            neg_map[idx_0] = weight;
            w0 += weight * weight;
        }
    }

    w0 = std::sqrt(w0);
    w1 = std::sqrt(w1);

    if (w0 > 1e-10)
        for (auto& [idx, w] : neg_map)
            w /= w0;
    if (w1 > 1e-10)
        for (auto& [idx, w] : pos_map)
            w /= w1;

    return {QRState(neg_map, state.n_bits), QRState(pos_map, state.n_bits), w0, w1};
}

void to_controlled_gate(const std::shared_ptr<QGate>& gate, uint32_t ctrl, bool phase,
                        std::vector<std::shared_ptr<QGate>>& out) {
    if (auto ry = std::dynamic_pointer_cast<RY>(gate)) {
        out.push_back(std::make_shared<CRY>(ctrl, phase, ry->theta, ry->target));
    } else if (auto cry = std::dynamic_pointer_cast<CRY>(gate)) {
        std::vector<uint32_t> ctrls  = {ctrl, cry->ctrl};
        std::vector<bool>     phases = {phase, cry->phase};
        out.push_back(std::make_shared<MCRY>(ctrls, phases, cry->theta, cry->target));
    } else if (auto cx = std::dynamic_pointer_cast<CX>(gate)) {
        std::vector<uint32_t> ctrls  = {ctrl, cx->ctrl};
        std::vector<bool>     phases = {phase, cx->phase};
        out.push_back(std::make_shared<MCRY>(ctrls, phases, M_PI, cx->target));
    } else if (auto mcry = std::dynamic_pointer_cast<MCRY>(gate)) {
        std::vector<uint32_t> ctrls  = {ctrl};
        std::vector<bool>     phases = {phase};
        ctrls.insert(ctrls.end(), mcry->ctrls.begin(), mcry->ctrls.end());
        phases.insert(phases.end(), mcry->phases.begin(), mcry->phases.end());
        out.push_back(std::make_shared<MCRY>(ctrls, phases, mcry->theta, mcry->target));
    } else if (auto x = std::dynamic_pointer_cast<X>(gate)) {
        out.push_back(std::make_shared<CX>(ctrl, phase, x->target));
    }
}

void prepare_dense_rec(const QRState& state, std::vector<std::shared_ptr<QGate>>& gates) {
    std::vector<uint32_t> supports;
    for (const auto& [index, weight] : state.index_to_weight)
        for (uint32_t i = 0; i < state.n_bits; i++)
            if (index & (1 << i)) {
                if (std::find(supports.begin(), supports.end(), i) == supports.end())
                    supports.push_back(i);
            }

    if (supports.size() <= 4 || state.cardinality() <= 100) {
        QCircuit circ = prepare_state(state, false);
        for (const auto& g : circ.pGates)
            gates.push_back(g);
        return;
    }

    uint32_t  pivot = select_pivot_qubit(state, supports);
    Cofactors cf    = compute_cofactors(state, pivot);

    gates.push_back(std::make_shared<RY>(pivot, 2.0 * std::atan2(cf.weight1, cf.weight0)));

    std::vector<std::shared_ptr<QGate>> pos_gates, neg_gates;
    if (cf.pos_state.cardinality() > 0)
        prepare_dense_rec(cf.pos_state, pos_gates);
    if (cf.neg_state.cardinality() > 0)
        prepare_dense_rec(cf.neg_state, neg_gates);

    for (const auto& g : pos_gates)
        to_controlled_gate(g, pivot, true, gates);
    for (const auto& g : neg_gates)
        to_controlled_gate(g, pivot, false, gates);
}

} // namespace

QCircuit prepare_state_dense(const QRState& state) {
    QCircuit                            circuit(state.n_bits);
    std::vector<std::shared_ptr<QGate>> gates;
    prepare_dense_rec(state, gates);
    for (auto it = gates.rbegin(); it != gates.rend(); ++it)
        circuit.add_gate(*it);
    return circuit;
}

} // namespace xyz
