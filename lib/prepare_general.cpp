#include "qcircuit.hpp"

#include <iostream>
#include <memory>
#include <optional>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xyz {
namespace {
struct bfs_state {
    uint32_t mem_idx;
    uint32_t cnot_cost;
    bool     operator<(const bfs_state& other) const { return cnot_cost > other.cnot_cost; }
    bfs_state(uint32_t mem_idx, uint32_t cnot_cost) : mem_idx(mem_idx), cnot_cost(cnot_cost) {}
};

struct memorized_state {
    uint32_t               prev;
    QRState                state;
    std::shared_ptr<QGate> gate;
    memorized_state(const QRState& state, uint32_t prev) : prev(prev), state(state) {}
    memorized_state(const QRState& state, uint32_t prev, std::shared_ptr<QGate> gate)
        : prev(prev), state(state), gate(gate) {}
};

std::vector<std::shared_ptr<QGate>> enumerate_gates(const QRState& state) {
    std::vector<std::shared_ptr<QGate>> gates;
    if (state.index_to_weight.size() == 1) {
        auto index = state.index_to_weight.begin()->first;
        for (uint32_t target = 0; target < state.n_bits; target++)
            if ((index >> target) & 1) {
                gates.push_back(std::make_shared<X>(target));
                return gates;
            }
    }
    for (uint32_t target = 0; target < state.n_bits; target++) {
        const auto            ry_table = state.to_ry_table(target);
        std::optional<double> theta;
        for (const auto& [index, t] : ry_table) {
            if (theta.has_value() && theta.value() != t) {
                theta.reset();
                break;
            }
            theta = t;
        }
        if (theta.has_value() && !Rotation::is_trivial(theta.value(), true)) {
            gates.push_back(std::make_shared<RY>(target, theta.value()));
            return gates;
        }
    }
    for (uint32_t target = 0; target < state.n_bits; target++) {
        const auto ry_table = state.to_ry_table(target);
        for (uint32_t ctrl = 0; ctrl < state.n_bits; ctrl++)
            for (bool phase : {false, true}) {
                if (ctrl == target)
                    continue;
                std::optional<double> theta;
                for (const auto& [index, t] : ry_table)
                    if ((bool)((index >> ctrl) & 1) == phase) {
                        if (theta.has_value() && theta.value() != t) {
                            theta.reset();
                            break;
                        }
                        theta = t;
                    }
                if (theta.has_value() && !Rotation::is_trivial(theta.value(), true)) {
                    gates.push_back(std::make_shared<CRY>(ctrl, phase, theta.value(), target));
                    gates.push_back(std::make_shared<CRY>(ctrl, phase, -M_PI + theta.value(), target));
                }
            }
    }
    for (uint32_t target = 0; target < state.n_bits; target++)
        for (uint32_t ctrl = 0; ctrl < state.n_bits; ctrl++) {
            if (ctrl == target)
                continue;
            gates.push_back(std::make_shared<CX>(ctrl, true, target));
        }
    return gates;
}

template <bool verbose> void prepare_state_impl(const QRState& state, QCircuit& circuit) {
    std::priority_queue<bfs_state>         q;
    std::vector<memorized_state>           states;
    std::unordered_map<uint64_t, uint32_t> cost;
    std::unordered_set<uint64_t>           visited;
    std::optional<uint32_t>                solution;

    states.push_back(memorized_state(state, (uint32_t)-1));
    q.push(bfs_state(0, 0));
    cost[state.repr()] = 0;

    while (!q.empty()) {
        auto e = q.top();
        q.pop();
        auto cur = states[e.mem_idx];
        auto r   = cur.state.repr();
        if (visited.find(r) != visited.end())
            continue;
        visited.insert(r);
        if (cur.state.is_ground()) {
            solution = e.mem_idx;
            break;
        }
        if constexpr (verbose)
            std::cout << "current_state: " << cur.state.to_string() << " repr = " << r << "\n";
        for (const auto& gate : enumerate_gates(cur.state)) {
            if constexpr (verbose)
                std::cout << "gate: " << *gate << "\n";
            auto new_state = (*gate)(cur.state, true);
            auto new_cost  = e.cnot_cost + gate->num_cnots();
            auto new_repr  = new_state.repr();
            if (visited.find(new_repr) != visited.end())
                continue;
            auto it = cost.find(new_repr);
            if (it == cost.end() || it->second > e.cnot_cost + 1) {
                cost[new_repr] = new_cost;
                q.push(bfs_state((uint32_t)states.size(), new_cost));
                states.push_back(memorized_state(new_state, e.mem_idx, gate));
            }
        }
    }

    if (!solution.has_value())
        return;

    for (uint32_t idx = solution.value(); idx != 0; idx = states[idx].prev)
        circuit.add_gate(states[idx].gate);
}
} // namespace

QCircuit prepare_state(const QRState& state, bool verbose) {
    QCircuit circuit(state.n_bits);
    if (verbose)
        prepare_state_impl<true>(state, circuit);
    else
        prepare_state_impl<false>(state, circuit);
    return circuit;
}

} // namespace xyz
