#pragma once

#include <algorithm>
#include <catch.hpp>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <unordered_set>
#include <vector>
#include "prepare-state.hpp"
#include "qcircuit.hpp"
#include "qgate.hpp"
#include "qstate.hpp"

namespace xyz::testutil {
inline QRState normalize(QRState s) {
    double ss = 0.0;
    for (const auto& [_, w] : s.index_to_weight)
        ss += w * w;
    REQUIRE(ss > 0.0);
    double inv = 1.0 / std::sqrt(ss);
    for (auto& [_, w] : s.index_to_weight)
        w *= inv;
    return s;
}

inline QRState make_state(uint32_t n, const std::vector<uint32_t>& indices, const std::vector<double>& weights) {
    REQUIRE(indices.size() == weights.size());
    std::map<uint32_t, double> m;
    for (size_t i = 0; i < indices.size(); i++)
        if (std::abs(weights[i]) >= QRState::eps)
            m[indices[i]] = weights[i];
    QRState s(m, n);
    return normalize(s);
}

inline QRState random_signed_sparse_state(uint32_t n, std::mt19937_64& rng, uint32_t max_support) {
    uint32_t                                dim = 1u << n;
    std::uniform_int_distribution<uint32_t> kdist(1u, std::min<uint32_t>(max_support, dim));
    std::uniform_int_distribution<uint32_t> idist(0u, dim - 1u);
    std::normal_distribution<double>        wdist(0.0, 1.0);
    uint32_t                                k = kdist(rng);
    std::unordered_set<uint32_t>            seen;
    std::vector<uint32_t>                   idx;
    std::vector<double>                     w;
    idx.reserve(k);
    w.reserve(k);
    while (idx.size() < k) {
        uint32_t x = idist(rng);
        if (seen.insert(x).second) {
            idx.push_back(x);
            w.push_back(wdist(rng));
        }
    }
    return make_state(n, idx, w);
}

inline void require_close(const QRState& a, const QRState& b, double eps = 1e-6) {
    REQUIRE(a.n_bits == b.n_bits);
    std::vector<uint32_t> keys;
    keys.reserve(a.index_to_weight.size() + b.index_to_weight.size());
    for (const auto& [k, _] : a.index_to_weight)
        keys.push_back(k);
    for (const auto& [k, _] : b.index_to_weight)
        keys.push_back(k);
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    double s = 1.0;
    for (auto k : keys) {
        double wa = 0.0, wb = 0.0;
        auto   ia = a.index_to_weight.find(k);
        auto   ib = b.index_to_weight.find(k);
        if (ia != a.index_to_weight.end())
            wa = ia->second;
        if (ib != b.index_to_weight.end())
            wb = ib->second;
        if (std::abs(wa) > eps && std::abs(wb) > eps) {
            if (wa * wb < 0.0)
                s = -1.0;
            break;
        }
    }
    for (auto k : keys) {
        double wa = 0.0, wb = 0.0;
        auto   ia = a.index_to_weight.find(k);
        auto   ib = b.index_to_weight.find(k);
        if (ia != a.index_to_weight.end())
            wa = ia->second;
        if (ib != b.index_to_weight.end())
            wb = ib->second;
        REQUIRE(std::abs(wa - s * wb) <= eps);
    }
}
} // namespace xyz::testutil
