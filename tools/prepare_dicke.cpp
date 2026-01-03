#include <cmdline.hpp>
#include <cstdint>
#include <iostream>
#include <map>
#include <xyz.hpp>

using namespace xyz;
using namespace cmdline;

parser CommandLineParser() {
    parser opt;
    opt.add<int>("n", 'n', "number of qubits", true);
    opt.add<int>("k", 'k', "number of excitations", true);
    return opt;
}

int main(int argc, char** argv) {
    auto opt = CommandLineParser();
    opt.parse_check(argc, argv);
    auto n       = opt.get<int>("n");
    auto k       = opt.get<int>("k");
    auto circuit = prepare_dicke_state(n, k);
    write_qasm2(circuit, "dicke_" + std::to_string(n) + "_" + std::to_string(k) + ".qasm");
    (void)simulate_circuit(circuit, ground_rstate(n), true);
    return 0;
}
