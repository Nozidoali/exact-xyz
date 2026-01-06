#include <cmdline.hpp>
#include <qcircuit.hpp>
#include <qgate.hpp>
#include <qstate.hpp>
#include <transpile.hpp>

using namespace xyz;
using namespace cmdline;

parser CommandLineParser() {
    parser opt;
    opt.add<int>("n", 'n', "number of qubits", true);
    opt.add<int>("k", 'k', "number of excitations", true);
    opt.add<std::string>("output", 'o', "output file", false, "");
    opt.add<double>("eps", 'e', "transpile epsilon", false, 1e-1);
    opt.add("transpile", 't', "transpile to Clifford+T");
    return opt;
}

int main(int argc, char** argv) {
    auto opt = CommandLineParser();
    opt.parse_check(argc, argv);
    auto n       = opt.get<int>("n");
    auto k       = opt.get<int>("k");
    auto circuit = prepare_dicke_state(n, k);

    if (opt.exist("transpile")) {
        double eps = opt.get<double>("eps");
        circuit    = transpile_clifford_t(circuit, eps);
    }

    if (opt.exist("output")) {
        std::string output = opt.get<std::string>("output");
        write_qasm2(circuit, output);
    } else {
        std::cout << circuit.to_qasm2();
    }

    std::cout << "CNOT count: " << circuit.num_cnots() << std::endl;

    return 0;
}
