#include <cstdint>
#include <iostream>
#include <map>
#include <xyz.hpp>

using namespace xyz;
int main() {
    uint32_t n_qubits  = 4;
    QRState  state_exp = dicke_state(n_qubits, 2);
    QCircuit qc        = prepare_sparse_state(state_exp);
    // qc = decompose_circuit( qc );
    write_qasm2(qc, "dicke_3_2.qasm");
    QRState state = simulate_circuit(qc, ground_rstate(n_qubits), true);
    return 0;
}