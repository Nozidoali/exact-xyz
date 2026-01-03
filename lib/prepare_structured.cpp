#include "qcircuit.hpp"

#include <array>
#include <cmath>
#include <queue>

namespace xyz {
QCircuit prepare_ghz(uint32_t n, bool log_depth) {
    QCircuit circuit(n);
    circuit.add_gate(std::make_shared<H>(0));
    if (log_depth)
        for (uint32_t i = 1u; i < n; i <<= 1)
            for (uint32_t j = 0; j < i && j + i < n; j++)
                circuit.add_gate(std::make_shared<CX>(j, true, j + i));
    else
        for (uint32_t i = 1; i < n; i++)
            circuit.add_gate(std::make_shared<CX>(0, true, i));
    return circuit;
}

QCircuit prepare_w(uint32_t n, bool log_depth, bool cnot_opt) {
    QCircuit circuit(n);
    circuit.add_gate(std::make_shared<X>(0));
    if (!log_depth) {
        for (uint32_t i = 1; i < n; i++) {
            uint32_t j     = i - 1;
            double   p     = 1.0 / (double)(n - j);
            double   theta = 2 * std::atan2(std::sqrt(1 - p), std::sqrt(p));
            if (cnot_opt) {
                circuit.add_gate(std::make_shared<RY>(i, -(theta - M_PI) / 2));
                circuit.add_gate(std::make_shared<CX>(j, true, i));
                circuit.add_gate(std::make_shared<RY>(i, (theta - M_PI) / 2));
            } else {
                circuit.add_gate(std::make_shared<CRY>(j, true, theta, i));
            }
            circuit.add_gate(std::make_shared<CX>(i, true, j));
        }
        return circuit;
    }

    std::queue<std::array<uint32_t, 3>> dicotomies;
    dicotomies.push({{0, n, n >> 1}});
    uint32_t q_next = 1;
    while (!dicotomies.empty()) {
        auto [q, total, curr] = dicotomies.front();
        dicotomies.pop();
        if (total < 2)
            continue;
        uint32_t total_l = total >> 1;
        uint32_t curr_l  = curr >> 1;
        uint32_t total_r = total - total_l;
        uint32_t curr_r  = curr - curr_l;
        if (total_l == 1 && curr_l == 1)
            dicotomies.push({q, total_r, curr_r});
        else {
            dicotomies.push({q, total_l, curr_l});
            dicotomies.push({q_next, total_r, curr_r});
        }
        double p     = (double)curr / (double)total;
        double theta = 2 * std::atan2(std::sqrt(1 - p), std::sqrt(p));
        if (cnot_opt) {
            circuit.add_gate(std::make_shared<RY>(q_next, -(theta - M_PI) / 2));
            circuit.add_gate(std::make_shared<CX>(q, true, q_next));
            circuit.add_gate(std::make_shared<RY>(q_next, (theta - M_PI) / 2));
        } else {
            circuit.add_gate(std::make_shared<CRY>(q, true, theta, q_next));
        }
        circuit.add_gate(std::make_shared<CX>(q_next, true, q));
        q_next++;
    }
    return circuit;
}
} // namespace xyz
