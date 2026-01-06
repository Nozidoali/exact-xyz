#include <cmdline.hpp>
#include <transpile.hpp>
#include <qcircuit.hpp>
#include <qgate.hpp>
#include <qstate.hpp>

using cmdline::parser;
using namespace xyz;

parser CommandLineParser() {
    parser opt;
    opt.add<std::string>("input", 'i', "path to input QASM2", false, "../data/input.qasm");
    opt.add<std::string>("output", 'o', "path to output QASM2", false, "out_ct.qasm");
    opt.add<double>("eps", 'e', "approximation tolerance", false, 1e-3);
    return opt;
}

int main(int argc, char** argv) {
    auto opt = CommandLineParser();
    opt.parse_check(argc, argv);

    auto in  = read_qasm2(opt.get<std::string>("input"));
    auto out = transpile_clifford_t(in, opt.get<double>("eps"));
    write_qasm2(out, opt.get<std::string>("output"));
    return 0;
}
