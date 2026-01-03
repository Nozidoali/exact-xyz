#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("prepare_sparse_state produces correct state", "[xyz]") {
    std::mt19937_64 rng(1);
    for (uint32_t n = 2; n <= 6; n++) {
        for (int it = 0; it < 50; it++) {
            auto    target = random_signed_sparse_state(n, rng, 6);
            auto    c      = prepare_sparse_state(target);
            QRState got    = simulate_circuit(c, ground_rstate(n), false);
            require_close(target, got);
        }
    }
}
