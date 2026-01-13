#pragma once

#include "qcircuit.hpp"

#include <memory>
#include <vector>

namespace xyz {

struct ReductionResult {
    QRState                             state;
    std::vector<std::shared_ptr<QGate>> gates;
};

ReductionResult support_reduction(const QRState& input_state);
ReductionResult cardinality_reduction_by_one(const QRState& state);

QCircuit prepare_state_auto(const QRState& state, bool verbose = false);
QCircuit prepare_state_dense(const QRState& state);
QCircuit prepare_sparse_state(const QRState& state);

} // namespace xyz
