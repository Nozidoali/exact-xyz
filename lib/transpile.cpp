#include "transpile.hpp"

#include "unitary.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace xyz {

QCircuit transpile_clifford_t(const QCircuit& in, double eps) {
    QCircuit lowered = decompose_circuit(in);
    QCircuit out(lowered.num_qbits);

    for (const auto& g : lowered.pGates) {
        if (auto ry = std::dynamic_pointer_cast<RY>(g)) {
            auto ry_gates = approximate_ry_sk(ry->target, ry->theta, eps);
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
