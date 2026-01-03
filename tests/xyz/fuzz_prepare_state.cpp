#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("prepare_state fuzz (random signed real)", "[xyz]") {
    std::mt19937_64 rng(2);
    for (uint32_t n = 2; n <= 3; n++) {
        for (int it = 0; it < 25; it++) {
            auto    target = random_signed_sparse_state(n, rng, 3);
            auto    c      = prepare_state(target, false);
            QRState got    = simulate_circuit(c, ground_rstate(n), false);
            require_close(target, got);
        }
    }
}
