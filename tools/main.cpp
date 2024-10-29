/*
 * Author: Hanyu Wang
 * Created time: 2024-03-30 18:35:38
 * Last Modified by: Hanyu Wang
 * Last Modified time: 2024-04-02 12:30:34
 */

#include <iostream>
#include <map>
#include <cstdint>
#include <qcircuit.hpp>

using namespace xyz;
int main() {
    QState state = dicke_state(5, 2);
    QCircuit qcircuit = prepare_state(state);
    std::cout << to_qasm2(qcircuit) << std::endl;
    return 0;
}