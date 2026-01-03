#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("prepare_ghz produces correct state", "[xyz]") {
    for (uint32_t n = 2; n <= 6; n++) {
        auto                       c   = prepare_ghz(n, false);
        QRState                    got = simulate_circuit(c, ground_rstate(n), false);
        double                     a   = 1.0 / std::sqrt(2.0);
        std::map<uint32_t, double> m;
        m[0]              = a;
        m[(1u << n) - 1u] = a;
        QRState exp(m, n);
        require_close(normalize(exp), got);
    }
}
