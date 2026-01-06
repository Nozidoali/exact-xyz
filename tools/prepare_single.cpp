#include <cmath>
#include <cmdline.hpp>
#include <iostream>
#include <qcircuit.hpp>
#include <qgate.hpp>
#include <qstate.hpp>
#include <transpile.hpp>

using namespace xyz;
using cmdline::parser;

parser CommandLineParser() {
    parser opt;
    opt.add<double>("theta", 't', "rotation angle for RY gate", false, M_PI / 4.0);
    opt.add<double>("eps", 'e', "transpile epsilon for Clifford+T", false, 1e-3);
    opt.add<std::string>("output", 'o', "output QASM file", false, "");
    return opt;
}

int main(int argc, char** argv) {
    auto opt = CommandLineParser();
    opt.parse_check(argc, argv);

    double theta = opt.get<double>("theta");
    double eps   = opt.get<double>("eps");

    // Create single-qubit circuit with RY rotation
    QCircuit circuit(1);
    circuit.add_gate(std::make_shared<RY>(0, theta));

    std::cout << "Original circuit with RY(" << theta << "):\n";
    std::cout << circuit.to_qasm2() << std::endl;

    // Transpile to Clifford+T
    circuit = transpile_clifford_t(circuit, eps);
    std::cout << "Clifford+T circuit (eps=" << eps << "):\n";
    std::cout << circuit.to_qasm2() << std::endl;

    if (opt.exist("output")) {
        std::string output = opt.get<std::string>("output");
        write_qasm2(circuit, output);
        std::cout << "Circuit written to: " << output << std::endl;
    }

    return 0;
}
