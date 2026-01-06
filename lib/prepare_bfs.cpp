#include "qcircuit.hpp"

#include <algorithm>
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
    uint32_t depth;
    bool     operator<(const bfs_state& other) const { return cnot_cost > other.cnot_cost; }
    bfs_state(uint32_t mem_idx, uint32_t cnot_cost, uint32_t depth)
        : mem_idx(mem_idx), cnot_cost(cnot_cost), depth(depth) {}
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

bool prepare_state_impl(const QRState& state, QCircuit& circuit, const bfs_params& params, bool verbose) {
    std::priority_queue<bfs_state>                     q;
    std::vector<memorized_state>                       states;
    std::unordered_map<QRState, uint32_t, QRStateHash> cost;
    std::unordered_set<QRState, QRStateHash>           visited;
    std::optional<uint32_t>                            solution;

    states.push_back(memorized_state(state, (uint32_t)-1));
    q.push(bfs_state(0, 0, 0));
    cost[state] = 0;

    while (!q.empty()) {
        auto e = q.top();
        q.pop();
        auto cur = states[e.mem_idx];
        if (visited.find(cur.state) != visited.end())
            continue;
        visited.insert(cur.state);
        if (cur.state.is_ground()) {
            solution = e.mem_idx;
            break;
        }
        if (e.depth >= params.max_depth)
            continue;
        if (verbose)
            std::cout << "current_state: " << cur.state.to_string() << "\n";
        auto gates = enumerate_gates(cur.state);
        if (gates.size() > params.max_neighbors) {
            std::sort(gates.begin(), gates.end(),
                      [](const auto& a, const auto& b) { return a->num_cnots() < b->num_cnots(); });
            gates.resize(params.max_neighbors);
        }
        for (const auto& gate : gates) {
            if (verbose)
                std::cout << "gate: " << *gate << "\n";
            auto new_state = (*gate)(cur.state, true);
            auto new_cost  = e.cnot_cost + gate->num_cnots();
            if (visited.find(new_state) != visited.end())
                continue;
            auto it = cost.find(new_state);
            if (it == cost.end() || it->second > e.cnot_cost + 1) {
                cost[new_state] = new_cost;
                q.push(bfs_state((uint32_t)states.size(), new_cost, e.depth + 1));
                states.push_back(memorized_state(new_state, e.mem_idx, gate));
            }
        }
    }

    if (!solution.has_value())
        return false;

    for (uint32_t idx = solution.value(); idx != 0; idx = states[idx].prev)
        circuit.add_gate(states[idx].gate);
    return true;
}
} // namespace

bool prepare_state_bfs(const QRState& state, QCircuit& circuit, const bfs_params& params, bool verbose) {
    return prepare_state_impl(state, circuit, params, verbose);
}

QCircuit prepare_state(const QRState& state, bool verbose) {
    QCircuit   circuit(state.n_bits);
    bfs_params params;
    prepare_state_impl(state, circuit, params, verbose);
    return circuit;
}

} // namespace xyz
