#include "qcircuit.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <optional>
#include <queue>
#include <vector>

namespace xyz {
namespace {

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
