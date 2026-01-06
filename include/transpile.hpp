#pragma once

#include "qcircuit.hpp"

#include <cstdint>
#include <vector>

namespace xyz {

QCircuit transpile_clifford_t(const QCircuit& in, double eps);

} // namespace xyz
