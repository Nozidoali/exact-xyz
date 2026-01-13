#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("prepare_state_dense simple correctness test", "[xyz][dense]") {
    SECTION("Bell state") {
        std::map<uint32_t, double> m = {{0, 1.0 / std::sqrt(2.0)}, {3, 1.0 / std::sqrt(2.0)}};
        QRState                    target(m, 2);
        auto                       c   = prepare_state_dense(target);
        QRState                    got = simulate_circuit(c, ground_rstate(2), false);
        require_close(target, got, 1e-3);
    }
    
    SECTION("W state") {
        double                     w = 1.0 / std::sqrt(3.0);
        std::map<uint32_t, double> m = {{1, w}, {2, w}, {4, w}};
        QRState                    target(m, 3);
        auto                       c   = prepare_state_dense(target);
        QRState                    got = simulate_circuit(c, ground_rstate(3), false);
        require_close(target, got, 1e-3);
    }
}

TEST_CASE("prepare_state_dense random states with relaxed tolerance", "[xyz][dense]") {
    SECTION("Small states with 0.1% tolerance") {
        for (uint32_t n = 3; n <= 5; n++) {
            for (int it = 0; it < 5; it++) {
                uint64_t seed  = 42 + n * 10 + it;
                auto     target = random_rstate(n, 4, seed);
                auto     c      = prepare_state_dense(target);
                QRState  got    = simulate_circuit(c, ground_rstate(n), false);
                require_close(target, got, 1e-3);
            }
        }
    }
}
