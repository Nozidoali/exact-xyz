#pragma once

#include "qgate.hpp"
#include "qstate.hpp"

#include <memory>
#include <vector>

namespace xyz {
struct bfs_params {
    uint32_t max_depth     = 12;
    uint32_t max_neighbors = 100;
    bfs_params()           = default;
    bfs_params(uint32_t max_depth, uint32_t max_neighbors) : max_depth(max_depth), max_neighbors(max_neighbors) {}
};

class QCircuit {
  public:
    uint32_t                            num_qbits = 0;
    std::vector<std::shared_ptr<QGate>> pGates;

  public:
    QCircuit() = default;
    QCircuit(uint32_t num_qbits) : num_qbits(num_qbits) {};
    void        add_gate(std::shared_ptr<QGate> gate);
    void        reverse();
    uint32_t    num_cnots() const;
    uint32_t    lev_cnots() const;
    std::string to_qasm2() const;
};

QCircuit decompose_circuit(const QCircuit& circuit);

QCircuit prepare_state(const QRState& state, bool verbose = false);
bool     prepare_state_bfs(const QRState& state, QCircuit& circuit, const bfs_params& params, bool verbose = false);
QCircuit prepare_ghz(uint32_t n, bool log_depth = false);
QCircuit prepare_w(uint32_t n, bool log_depth = false, bool cnot_opt = false);
QCircuit prepare_dicke_state(int n, int k);

QCircuit resyn(const QCircuit& circuit, bool verbose = false);

void     write_qasm2(const QCircuit& circuit, const std::string& filename);
QCircuit read_qasm2(const std::string& filename, bool verbose = false);

QRState simulate_circuit(const QCircuit& circuit, const QRState& state, bool verbose = false);

QCircuit transpile_clifford_t(const QCircuit& in, double eps);

} // namespace xyz
