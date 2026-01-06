#include <cmdline.hpp>
#include <cstdint>
#include <transpile.hpp>
#include <xyz.hpp>

using cmdline::parser;
using namespace xyz;

namespace {
struct Counts {
    uint64_t cx  = 0;
    uint64_t t   = 0;
    uint64_t tdg = 0;
    uint64_t s   = 0;
    uint64_t sdg = 0;
    uint64_t x   = 0;
    uint64_t z   = 0;
};

Counts count_gates(const QCircuit& c) {
    Counts out;
    for (const auto& g : c.pGates) {
        if (std::dynamic_pointer_cast<CX>(g))
            out.cx++;
        else if (std::dynamic_pointer_cast<Tdg>(g))
            out.tdg++;
        else if (std::dynamic_pointer_cast<T>(g))
            out.t++;
        else if (std::dynamic_pointer_cast<Sdg>(g))
            out.sdg++;
        else if (std::dynamic_pointer_cast<S>(g))
            out.s++;
        else if (std::dynamic_pointer_cast<X>(g))
            out.x++;
        else if (std::dynamic_pointer_cast<Z>(g))
            out.z++;
    }
    return out;
}
} // namespace

parser CommandLineParser() {
    parser opt;
    opt.add<int>("n", 'n', "number of qubits", true);
    opt.add<int>("cardinality", 'c', "target state support size", false, 8);
    opt.add<uint64_t>("seed", 's', "seed (0 = nondeterministic)", false, 0);
    opt.add<double>("eps", 'e', "transpile epsilon", false, 1e-3);
    opt.add<std::string>("out", 'o', "output prefix (writes *_prep.qasm and *_ct.qasm)", false, "");
    opt.add("no_transpile", 0, "skip Clifford+T transpilation");
    opt.add("json", 0, "print JSON");
    opt.add("verbose", 'v', "verbose output");
    return opt;
}

int main(int argc, char** argv) {
    auto opt = CommandLineParser();
    opt.parse_check(argc, argv);

    uint32_t n = (uint32_t)opt.get<int>("n");
    uint32_t c = (uint32_t)opt.get<int>("cardinality");
    uint64_t s = opt.get<uint64_t>("seed");
    double   e = opt.get<double>("eps");
    bool     v = opt.exist("verbose");

    auto target = random_rstate(n, c, s);
    auto prep   = prepare_state_auto(target, v);
    auto prep_d = decompose_circuit(prep);

    bool     do_transpile = !opt.exist("no_transpile");
    QCircuit ct;
    if (do_transpile)
        ct = transpile_clifford_t(prep, e);

    auto   prep_counts = count_gates(prep_d);
    Counts ct_counts;
    if (do_transpile)
        ct_counts = count_gates(ct);

    double ct_err_max = 0.0;
    if (do_transpile) {
        auto lowered = decompose_circuit(prep);
        for (const auto& g : lowered.pGates) {
            if (auto ry = std::dynamic_pointer_cast<RY>(g)) {
                auto   word = sk::synthesize_rz(ry->theta, e);
                double err  = sk::dist(sk::word_matrix(word), sk::rz_matrix(ry->theta));
                if (err > ct_err_max)
                    ct_err_max = err;
            }
        }
    }

    auto out_prefix = opt.get<std::string>("out");
    if (!out_prefix.empty()) {
        write_qasm2(prep_d, out_prefix + "_prep.qasm");
        if (do_transpile)
            write_qasm2(ct, out_prefix + "_ct.qasm");
    }

    bool json = opt.exist("json");
    if (json) {
        std::cout << "{";
        std::cout << "\"n\":" << n << ",\"cardinality\":" << c << ",\"seed\":" << s << ",\"eps\":" << e;
        std::cout << ",\"prep\":{\"cx\":" << prep_counts.cx << ",\"t\":" << prep_counts.t
                  << ",\"tdg\":" << prep_counts.tdg << ",\"s\":" << prep_counts.s << ",\"sdg\":" << prep_counts.sdg
                  << ",\"x\":" << prep_counts.x << ",\"z\":" << prep_counts.z << "}";
        if (do_transpile) {
            std::cout << ",\"ct\":{\"cx\":" << ct_counts.cx << ",\"t\":" << ct_counts.t << ",\"tdg\":" << ct_counts.tdg
                      << ",\"s\":" << ct_counts.s << ",\"sdg\":" << ct_counts.sdg << ",\"x\":" << ct_counts.x
                      << ",\"z\":" << ct_counts.z << "}";
            std::cout << ",\"ct_err_max\":" << ct_err_max;
        }
        std::cout << "}\n";
    } else {
        std::cout << "prep: cx=" << prep_counts.cx << " t=" << prep_counts.t << " tdg=" << prep_counts.tdg
                  << " s=" << prep_counts.s << " sdg=" << prep_counts.sdg << " x=" << prep_counts.x
                  << " z=" << prep_counts.z << "\n";
        if (do_transpile)
            std::cout << "ct:   cx=" << ct_counts.cx << " t=" << ct_counts.t << " tdg=" << ct_counts.tdg
                      << " s=" << ct_counts.s << " sdg=" << ct_counts.sdg << " x=" << ct_counts.x
                      << " z=" << ct_counts.z << "\n";
    }

    return 0;
}
