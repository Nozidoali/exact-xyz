#include <cmdline.hpp>
#include <cstdint>
#include <qcircuit.hpp>
#include <qgate.hpp>
#include <qstate.hpp>

using namespace xyz;
using cmdline::parser;

parser CommandLineParser() {
    parser opt;
    opt.add<int>("n", 'n', "number of qubits", true);
    opt.add<int>("cardinality", 'c', "number of non-zero amplitudes in the target state", false, 8);
    opt.add<uint64_t>("seed", 's', "RNG seed (0 = non-deterministic)", false, 0);
    opt.add<std::string>("output", 'o', "output QASM2 file", false, "sparse_state.qasm");
    return opt;
}

int main(int argc, char** argv) {
    auto opt = CommandLineParser();
    opt.parse_check(argc, argv);

    const uint32_t n = (uint32_t)opt.get<int>("n");
    const uint32_t c = (uint32_t)opt.get<int>("cardinality");
    const uint64_t s = opt.get<uint64_t>("seed");

    QRState  target = random_rstate(n, c, s);
    QCircuit qc     = prepare_sparse_state(target);

    write_qasm2(qc, opt.get<std::string>("output"));
    (void)simulate_circuit(qc, ground_rstate(n), true);
    return 0;
}