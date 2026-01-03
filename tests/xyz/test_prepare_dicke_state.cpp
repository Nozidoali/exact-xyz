#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("prepare_dicke_state produces correct state", "[xyz]") {
    for (uint32_t n = 1; n <= 6; n++) {
        for (uint32_t k = 0; k <= n; k++) {
            auto    c   = prepare_dicke_state((int)n, (int)k);
            QRState got = simulate_circuit(c, ground_rstate(n), false);
            QRState exp = dicke_state(n, k);
            require_close(exp, got);
        }
    }
}
