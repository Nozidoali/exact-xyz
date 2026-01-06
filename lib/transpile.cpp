#include "transpile.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>

namespace xyz::sk {
namespace {

constexpr double kPi       = 3.141592653589793238462643383279502884;
constexpr double kSqrt2Inv = 0.707106781186547524400844362104849039;

SU2 from_array(const std::array<std::complex<double>, 4>& m) {
    return SU2{m[0], m[1], m[2], m[3]};
}

std::complex<double> det(const SU2& u) {
    return u.a * u.d - u.b * u.c;
}

SU2 scale(const SU2& u, std::complex<double> s) {
    return SU2{u.a * s, u.b * s, u.c * s, u.d * s};
}

struct Key {
    int32_t a_re, a_im, b_re, b_im;
    bool    operator==(const Key& o) const { return a_re == o.a_re && a_im == o.a_im && b_re == o.b_re && b_im == o.b_im; }
};

struct KeyHash {
    size_t operator()(const Key& k) const noexcept {
        uint64_t h   = 1469598103934665603ull;
        auto     mix = [&](uint64_t x) {
            h ^= x;
            h *= 1099511628211ull;
        };
        mix((uint32_t)k.a_re);
        mix((uint32_t)k.a_im);
        mix((uint32_t)k.b_re);
        mix((uint32_t)k.b_im);
        return (size_t)h;
    }
};

Key quantize(const SU2& u, double step) {
    auto q = [&](double x) -> int32_t { return (int32_t)std::llround(x / step); };
    return Key{q(std::real(u.a)), q(std::imag(u.a)), q(std::real(u.b)), q(std::imag(u.b))};
}

std::vector<Gate> build_left_word(uint32_t kL, uint32_t bitsL) {
    std::vector<Gate> w;
    w.reserve(2 * kL + 1);
    if (bitsL & 1u)
        w.push_back(Gate::H);
    for (uint32_t i = 0; i < kL; ++i) {
        w.push_back(Gate::T);
        if ((bitsL >> (i + 1)) & 1u)
            w.push_back(Gate::H);
    }
    return w;
}

std::vector<Gate> build_right_word(uint32_t kT, uint32_t kL, uint32_t bitsR) {
    std::vector<Gate> w;
    w.reserve(2 * (kT - kL));
    for (uint32_t i = kL; i < kT; ++i) {
        w.push_back(Gate::T);
        if ((bitsR >> (i - kL)) & 1u)
            w.push_back(Gate::H);
    }
    return w;
}

} // namespace

SU2 identity() {
    return SU2{1, 0, 0, 1};
}

SU2 mul(const SU2& x, const SU2& y) {
    return SU2{
        x.a * y.a + x.b * y.c,
        x.a * y.b + x.b * y.d,
        x.c * y.a + x.d * y.c,
        x.c * y.b + x.d * y.d,
    };
}

SU2 dagger(const SU2& u) {
    return SU2{std::conj(u.d), -std::conj(u.b), -std::conj(u.c), std::conj(u.a)};
}

SU2 normalize_det(const SU2& u) {
    auto d = det(u);
    if (std::abs(d) < 1e-14)
        return u;
    auto phase = std::sqrt(d);
    return scale(u, 1.0 / phase);
}

double dist(const SU2& u, const SU2& v) {
    auto   uv    = normalize_det(mul(dagger(u), v));
    double tr_re = std::real(uv.a + uv.d) * 0.5;
    tr_re        = std::clamp(tr_re, -1.0, 1.0);
    return std::acos(tr_re);
}

SU2 gate_matrix(Gate g) {
    using cd = std::complex<double>;
    using namespace std::complex_literals;
    if (g == Gate::H) {
        return normalize_det(from_array({cd(kSqrt2Inv, 0), cd(kSqrt2Inv, 0), cd(kSqrt2Inv, 0), cd(-kSqrt2Inv, 0)}));
    }
    if (g == Gate::T) {
        return normalize_det(from_array({cd(1, 0), cd(0, 0), cd(0, 0), std::exp(1i * kPi / 4.0)}));
    }
    return normalize_det(from_array({cd(1, 0), cd(0, 0), cd(0, 0), std::exp(-1i * kPi / 4.0)}));
}

SU2 rz_matrix(double theta) {
    using namespace std::complex_literals;
    return SU2{std::exp(-1i * theta / 2.0), 0, 0, std::exp(1i * theta / 2.0)};
}

SU2 word_matrix(const std::vector<Gate>& w) {
    SU2 u = identity();
    for (auto g : w)
        u = mul(u, gate_matrix(g));
    return u;
}

std::vector<Gate> invert_word(const std::vector<Gate>& w) {
    std::vector<Gate> out;
    out.reserve(w.size());
    for (auto it = w.rbegin(); it != w.rend(); ++it) {
        if (*it == Gate::T)
            out.push_back(Gate::Tdg);
        else if (*it == Gate::Tdg)
            out.push_back(Gate::T);
        else
            out.push_back(*it);
    }
    return out;
}

std::vector<Gate> synthesize_rz(double theta, double eps) {
    SU2               target = rz_matrix(theta);
    std::vector<Gate> best;
    double            best_d = std::numeric_limits<double>::infinity();

    const int k_max = (eps <= 1e-3) ? 28 : (eps <= 1e-2) ? 22 : 18;
    for (int kT = 0; kT <= k_max; ++kT) {
        int kL         = kT / 2;
        int left_bits  = kL + 1;
        int right_bits = kT - kL;

        uint32_t left_bits_n  = (uint32_t)(1u << left_bits);
        uint32_t right_bits_n = (right_bits == 0) ? 1u : (uint32_t)(1u << right_bits);

        double                                                                        step = 0.25;
        std::unordered_map<Key, std::vector<std::pair<SU2, std::vector<Gate>>>, KeyHash> buckets;
        buckets.reserve(right_bits_n * 2);

        for (uint32_t mask = 0; mask < right_bits_n; ++mask) {
            auto w = build_right_word((uint32_t)kT, (uint32_t)kL, mask);
            auto u = word_matrix(w);
            buckets[quantize(u, step)].push_back({u, std::move(w)});
        }

        for (uint32_t mask = 0; mask < left_bits_n; ++mask) {
            auto wL = build_left_word((uint32_t)kL, mask);
            auto uL = word_matrix(wL);

            auto v  = mul(target, dagger(uL));
            Key  kq = quantize(v, step);

            for (int da = -2; da <= 2; ++da)
                for (int db = -2; db <= 2; ++db)
                    for (int dc = -2; dc <= 2; ++dc)
                        for (int dd = -2; dd <= 2; ++dd) {
                            Key  kk{kq.a_re + da, kq.a_im + db, kq.b_re + dc, kq.b_im + dd};
                            auto it = buckets.find(kk);
                            if (it == buckets.end())
                                continue;
                            for (const auto& cand : it->second) {
                                const auto& uR = cand.first;
                                auto        u  = mul(uR, uL);
                                double      d  = dist(u, target);
                                if (d < best_d) {
                                    best_d = d;
                                    best   = wL;
                                    best.insert(best.end(), cand.second.begin(), cand.second.end());
                                }
                            }
                        }
        }

        if (best_d <= eps)
            break;
    }

    return best;
}

} // namespace xyz::sk

namespace xyz {
namespace {

void append_word_ht(QCircuit& out, uint32_t target, const std::vector<sk::Gate>& word) {
    for (auto g : word) {
        if (g == sk::Gate::H) {
            out.add_gate(std::make_shared<H>(target));
        } else if (g == sk::Gate::T) {
            out.add_gate(std::make_shared<T>(target));
        } else {
            out.add_gate(std::make_shared<Tdg>(target));
        }
    }
}

} // namespace

QCircuit transpile_clifford_t(const QCircuit& in, double eps) {
    QCircuit lowered = decompose_circuit(in);
    QCircuit out(lowered.num_qbits);

    for (const auto& g : lowered.pGates) {
        if (auto ry = std::dynamic_pointer_cast<RY>(g)) {
            auto rz_word = sk::synthesize_rz(ry->theta, eps);
            out.add_gate(std::make_shared<S>(ry->target));
            out.add_gate(std::make_shared<H>(ry->target));
            append_word_ht(out, ry->target, rz_word);
            out.add_gate(std::make_shared<H>(ry->target));
            out.add_gate(std::make_shared<Sdg>(ry->target));
            continue;
        }
        out.add_gate(g);
    }

    return out;
}

} // namespace xyz

