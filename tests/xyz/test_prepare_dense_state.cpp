#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("prepare_state_dense produces correct state for random states", "[xyz]") {
    auto make_positive_random_state = [](uint32_t n, uint32_t card, uint64_t seed) {
        auto   state = random_rstate(n, card, seed);
        double norm  = 0.0;
        for (auto& [idx, w] : state.index_to_weight) {
            w = std::abs(w);
            norm += w * w;
        }
        norm = std::sqrt(norm);
        for (auto& [idx, w] : state.index_to_weight) {
            w /= norm;
        }
        return state;
    };
    
    SECTION("Small states") {
        for (uint32_t n = 3; n <= 5; n++) {
            for (int it = 0; it < 20; it++) {
                auto    target = make_positive_random_state(n, 8, 42 + n * 100 + it);
                auto    c      = prepare_state_dense(target);
                QRState got    = simulate_circuit(c, ground_rstate(n), false);
                require_close(target, got, 1e-4);
            }
        }
    }
    
    SECTION("Medium states") {
        for (uint32_t n = 6; n <= 7; n++) {
            for (int it = 0; it < 10; it++) {
                auto    target = make_positive_random_state(n, 16, 1000 + n * 100 + it);
                auto    c      = prepare_state_dense(target);
                QRState got    = simulate_circuit(c, ground_rstate(n), false);
                require_close(target, got, 1e-4);
            }
        }
    }
    
    SECTION("Larger states") {
        for (int it = 0; it < 5; it++) {
            auto    target = make_positive_random_state(8, 32, 2000 + it);
            auto    c      = prepare_state_dense(target);
            QRState got    = simulate_circuit(c, ground_rstate(8), false);
            require_close(target, got, 1e-4);
        }
    }
}

TEST_CASE("prepare_state_dense with fixed random states", "[xyz]") {
    auto abs_normalize = [](const QRState& state) {
        std::map<uint32_t, double> m;
        double                     norm = 0.0;
        for (const auto& [idx, w] : state.index_to_weight) {
            double abs_w = std::abs(w);
            m[idx]       = abs_w;
            norm += abs_w * abs_w;
        }
        norm = std::sqrt(norm);
        for (auto& [idx, w] : m) {
            w /= norm;
        }
        return QRState(m, state.n_bits);
    };
    
    SECTION("n=4 cardinality=8") {
        auto    target = abs_normalize(random_rstate(4, 8, 42));
        auto    c      = prepare_state_dense(target);
        QRState got    = simulate_circuit(c, ground_rstate(4), false);
        require_close(target, got, 1e-4);
    }
    
    SECTION("n=5 cardinality=16") {
        auto    target = abs_normalize(random_rstate(5, 16, 100));
        auto    c      = prepare_state_dense(target);
        QRState got    = simulate_circuit(c, ground_rstate(5), false);
        require_close(target, got, 1e-4);
    }
    
    SECTION("n=6 cardinality=32") {
        auto    target = abs_normalize(random_rstate(6, 32, 200));
        auto    c      = prepare_state_dense(target);
        QRState got    = simulate_circuit(c, ground_rstate(6), false);
        require_close(target, got, 1e-4);
    }
}

TEST_CASE("prepare_state_dense edge cases", "[xyz]") {
    SECTION("Single amplitude") {
        std::map<uint32_t, double> m = {{5, 1.0}};
        QRState                    target(m, 4);
        auto                       c   = prepare_state_dense(target);
        QRState                    got = simulate_circuit(c, ground_rstate(4), false);
        require_close(target, got, 1e-4);
    }
    
    SECTION("Two amplitudes") {
        std::map<uint32_t, double> m = {{0, 0.6}, {15, 0.8}};
        QRState                    target(m, 4);
        auto                       c   = prepare_state_dense(target);
        QRState                    got = simulate_circuit(c, ground_rstate(4), false);
        require_close(target, got, 1e-4);
    }
    
    SECTION("Uniform superposition") {
        std::map<uint32_t, double> m;
        for (uint32_t i = 0; i < 8; i++) {
            m[i] = 1.0 / std::sqrt(8.0);
        }
        QRState target(m, 3);
        auto    c   = prepare_state_dense(target);
        QRState got = simulate_circuit(c, ground_rstate(3), false);
        require_close(target, got, 1e-4);
    }
    
    SECTION("Dense state all amplitudes") {
        std::mt19937_64                   rng(123);
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::map<uint32_t, double>        m;
        double                            norm = 0.0;
        for (uint32_t i = 0; i < 16; i++) {
            double w = std::abs(dist(rng));
            m[i]     = w;
            norm += w * w;
        }
        norm = std::sqrt(norm);
        for (auto& [_, w] : m) {
            w /= norm;
        }
        QRState target(m, 4);
        auto    c   = prepare_state_dense(target);
        QRState got = simulate_circuit(c, ground_rstate(4), false);
        require_close(target, got, 1e-4);
    }
}

TEST_CASE("prepare_state_dense performance test", "[xyz][.performance]") {
    SECTION("Varying cardinality") {
        std::vector<uint32_t> cardinalities = {4, 8, 16, 32, 50};
        for (size_t i = 0; i < cardinalities.size(); i++) {
            auto    target = random_rstate(6, cardinalities[i], 999 + i);
            auto    c      = prepare_state_dense(target);
            QRState got    = simulate_circuit(c, ground_rstate(6), false);
            INFO("n=6, card=" << cardinalities[i] << ", gates=" << c.pGates.size());
            require_close(target, got);
        }
    }
    
    SECTION("Varying qubits") {
        for (uint32_t n = 3; n <= 8; n++) {
            auto    target = random_rstate(n, 8, 5000 + n);
            auto    c      = prepare_state_dense(target);
            QRState got    = simulate_circuit(c, ground_rstate(n), false);
            INFO("n=" << n << ", card=8, gates=" << c.pGates.size());
            require_close(target, got);
        }
    }
}

TEST_CASE("prepare_state_dense with various seeds", "[xyz]") {
    auto abs_normalize = [](const QRState& s) {
        std::map<uint32_t, double> m;
        double                     norm = 0.0;
        for (const auto& [idx, w] : s.index_to_weight) {
            double abs_w = std::abs(w);
            m[idx]       = abs_w;
            norm += abs_w * abs_w;
        }
        norm = std::sqrt(norm);
        for (auto& [idx, w] : m) {
            w /= norm;
        }
        return QRState(m, s.n_bits);
    };
    
    SECTION("Multiple seeds") {
        std::vector<uint64_t> seeds = {1, 42, 100, 1000, 12345};
        for (auto seed : seeds) {
            auto    target = abs_normalize(random_rstate(5, 10, seed));
            auto    c      = prepare_state_dense(target);
            QRState got    = simulate_circuit(c, ground_rstate(5), false);
            require_close(target, got, 1e-4);
        }
    }
}
