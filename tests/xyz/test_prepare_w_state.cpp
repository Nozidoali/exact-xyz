#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

static QRState expected_w(uint32_t n) {
    std::map<uint32_t, double> m;
    double                     a = 1.0 / std::sqrt((double)n);
    for (uint32_t i = 0; i < n; i++)
        m[1u << i] = a;
    QRState s(m, n);
    return normalize(s);
}

TEST_CASE("prepare_w produces correct state", "[xyz]") {
    for (uint32_t n = 2; n <= 6; n++) {
        auto exp = expected_w(n);
        for (bool log_depth : {false, true}) {
            for (bool cnot_opt : {false, true}) {
                DYNAMIC_SECTION("n=" << n << " log=" << log_depth << " opt=" << cnot_opt) {
                    auto    c   = prepare_w(n, log_depth, cnot_opt);
                    QRState got = simulate_circuit(c, ground_rstate(n), false);
                    INFO(got.to_string());
                    require_close(exp, got);
                }
            }
        }
    }
}
