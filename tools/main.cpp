#include <iostream>
#include <map>
#include <cstdint>
#include <qcircuit.hpp>

using namespace xyz;
int main() {
    QCircuit qc = read_qasm2("output.qasm");
    std::cout << qc.num_cnots() << std::endl;
    QCircuit new_qc = resyn(qc);
    std::cout << new_qc.num_cnots() << std::endl;

    // QState state = dicke_state(5, 2);
    // QCircuit qcircuit = prepare_state(state);
    // qcircuit = decompose_circuit(qcircuit);
    // std::cout << qcircuit.num_cnots() << std::endl;
    // write_qasm2(qcircuit, "output.qasm");
    return 0;
}