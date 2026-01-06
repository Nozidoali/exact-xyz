#include "transpile.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace xyz {
namespace {

std::vector<std::shared_ptr<QGate>> approximate_ry(uint32_t target, double theta, double eps) {
    std::vector<std::shared_ptr<QGate>> gates;

    double base_step = M_PI / 4.0;
    int    max_steps = (int)(std::abs(theta) / base_step) + 10;

    double best_angle = 0.0;
    double best_error = std::abs(theta);
    int    best_n     = 0;

    for (int n = -max_steps; n <= max_steps; n++) {
        double angle = n * base_step;
        double error = std::abs(theta - angle);
        if (error < best_error) {
            best_error = error;
            best_angle = angle;
            best_n     = n;
        }
        if (error <= eps)
            break;
    }

    for (int i = 0; i < std::abs(best_n); i++) {
        gates.push_back(std::make_shared<H>(target));
        if (best_n > 0)
            gates.push_back(std::make_shared<T>(target));
        else
            gates.push_back(std::make_shared<Tdg>(target));
        gates.push_back(std::make_shared<H>(target));
    }

    return gates;
}

} // namespace

QCircuit transpile_clifford_t(const QCircuit& in, double eps) {
    QCircuit lowered = decompose_circuit(in);
    QCircuit out(lowered.num_qbits);

    for (const auto& g : lowered.pGates) {
        if (auto ry = std::dynamic_pointer_cast<RY>(g)) {
            auto ry_gates = approximate_ry(ry->target, ry->theta, eps);
            for (const auto& gate : ry_gates) {
                out.add_gate(gate);
            }
            continue;
        }
        out.add_gate(g);
    }

    return out;
}

} // namespace xyz
