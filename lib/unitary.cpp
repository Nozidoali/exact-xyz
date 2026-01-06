#include "unitary.hpp"

#include "qgate.hpp"

#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace xyz {

Matrix2x2 mmul(const Matrix2x2& a, const Matrix2x2& b) {
    Matrix2x2 r{};
    r[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
    r[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
    r[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
    r[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];
    return r;
}

Matrix2x2 dagger(const Matrix2x2& a) {
    Matrix2x2 r{};
    r[0][0] = std::conj(a[0][0]);
    r[0][1] = std::conj(a[1][0]);
    r[1][0] = std::conj(a[0][1]);
    r[1][1] = std::conj(a[1][1]);
    return r;
}

Complex tr(const Matrix2x2& a) {
    return a[0][0] + a[1][1];
}

Complex det2(const Matrix2x2& a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

Matrix2x2 normalize_to_su2(const Matrix2x2& u) {
    Complex d = det2(u);
    if (std::abs(d) < 1e-15)
        return u;
    Complex s = std::sqrt(d);
    if (std::abs(s) < 1e-15)
        return u;
    Matrix2x2 r = u;
    r[0][0] /= s;
    r[0][1] /= s;
    r[1][0] /= s;
    r[1][1] /= s;
    return r;
}

double dist_phase_invariant(const Matrix2x2& U, const Matrix2x2& V) {
    Matrix2x2 W = mmul(dagger(U), V);
    double    t = std::abs(tr(W)) * 0.5;
    if (t > 1.0)
        t = 1.0;
    if (t < 0.0)
        t = 0.0;
    return 1.0 - t;
}

Matrix2x2 mat_I() {
    return Matrix2x2{{{Complex(1, 0), Complex(0, 0)}, {Complex(0, 0), Complex(1, 0)}}};
}

Matrix2x2 mat_H() {
    double s = 1.0 / std::sqrt(2.0);
    return Matrix2x2{{{Complex(s, 0), Complex(s, 0)}, {Complex(s, 0), Complex(-s, 0)}}};
}

Matrix2x2 mat_T() {
    double a = M_PI / 4.0;
    return Matrix2x2{{{Complex(1, 0), Complex(0, 0)}, {Complex(0, 0), std::exp(Complex(0, a))}}};
}

Matrix2x2 mat_Tdg() {
    double a = -M_PI / 4.0;
    return Matrix2x2{{{Complex(1, 0), Complex(0, 0)}, {Complex(0, 0), std::exp(Complex(0, a))}}};
}

Matrix2x2 mat_Ry(double theta) {
    double c = std::cos(theta * 0.5);
    double s = std::sin(theta * 0.5);
    return Matrix2x2{{{Complex(c, 0), Complex(-s, 0)}, {Complex(s, 0), Complex(c, 0)}}};
}

namespace {

Matrix2x2 mat_X() {
    return Matrix2x2{{{Complex(0, 0), Complex(1, 0)}, {Complex(1, 0), Complex(0, 0)}}};
}
Matrix2x2 mat_Y() {
    return Matrix2x2{{{Complex(0, 0), Complex(0, -1)}, {Complex(0, 1), Complex(0, 0)}}};
}
Matrix2x2 mat_Z() {
    return Matrix2x2{{{Complex(1, 0), Complex(0, 0)}, {Complex(0, 0), Complex(-1, 0)}}};
}

Matrix2x2 exp_i_alpha_sigma(double nx, double ny, double nz, double alpha) {
    double    ca = std::cos(alpha);
    double    sa = std::sin(alpha);
    Matrix2x2 nS{};
    Matrix2x2 X = mat_X(), Y = mat_Y(), Z = mat_Z();
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) {
            nS[r][c] = nx * X[r][c] + ny * Y[r][c] + nz * Z[r][c];
        }
    Matrix2x2 I = mat_I();
    Matrix2x2 U{};
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) {
            U[r][c] = ca * I[r][c] + Complex(0, 1) * sa * nS[r][c];
        }
    return U;
}

void su2_axis_angle(const Matrix2x2& U_in, double& nx, double& ny, double& nz, double& phi) {
    Matrix2x2 U    = normalize_to_su2(U_in);
    double    cphi = std::real(tr(U)) * 0.5;
    if (cphi > 1.0)
        cphi = 1.0;
    if (cphi < -1.0)
        cphi = -1.0;
    phi         = std::acos(cphi);
    double sphi = std::sin(phi);
    if (sphi < 1e-12) {
        nx = 0;
        ny = 0;
        nz = 1;
        return;
    }
    Complex U12 = U[0][1];
    Complex U11 = U[0][0];
    ny          = std::real(U12) / sphi;
    nx          = std::imag(U12) / sphi;
    nz          = std::imag(U11) / sphi;

    double nrm = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (nrm < 1e-12) {
        nx = 0;
        ny = 0;
        nz = 1;
    } else {
        nx /= nrm;
        ny /= nrm;
        nz /= nrm;
    }
}

Matrix2x2 commutator(const Matrix2x2& V, const Matrix2x2& W) {
    Matrix2x2 Vinv = dagger(V);
    Matrix2x2 Winv = dagger(W);
    return mmul(mmul(mmul(V, W), Vinv), Winv);
}

Matrix2x2 align_z_to_n(double nx, double ny, double nz) {
    if (std::abs(nx) < 1e-12 && std::abs(ny) < 1e-12) {
        if (nz >= 0)
            return mat_I();
        return exp_i_alpha_sigma(1, 0, 0, M_PI / 2.0);
    }
    double ax = -ny, ay = nx, az = 0.0;
    double an = std::sqrt(ax * ax + ay * ay + az * az);
    ax /= an;
    ay /= an;
    az /= an;
    double gamma = std::acos(std::max(-1.0, std::min(1.0, nz)));
    return exp_i_alpha_sigma(ax, ay, az, gamma * 0.5);
}

void decompose_near_identity_commutator(const Matrix2x2& Delta_in, Matrix2x2& V, Matrix2x2& W) {
    double nx, ny, nz, phi;
    su2_axis_angle(Delta_in, nx, ny, nz, phi);

    if (phi < 1e-10) {
        V = mat_I();
        W = mat_I();
        return;
    }

    auto K_of = [](double alpha) {
        Matrix2x2 V0 = exp_i_alpha_sigma(1, 0, 0, alpha);
        Matrix2x2 W0 = exp_i_alpha_sigma(0, 1, 0, alpha);
        return commutator(V0, W0);
    };

    double lo = 0.0, hi = 0.6;
    for (int it = 0; it < 60; it++) {
        double mid = 0.5 * (lo + hi);
        double kx, ky, kz, kphi;
        su2_axis_angle(K_of(mid), kx, ky, kz, kphi);
        if (kphi < phi)
            lo = mid;
        else
            hi = mid;
    }
    double alpha = 0.5 * (lo + hi);

    Matrix2x2 V0 = exp_i_alpha_sigma(1, 0, 0, alpha);
    Matrix2x2 W0 = exp_i_alpha_sigma(0, 1, 0, alpha);

    Matrix2x2 A    = align_z_to_n(nx, ny, nz);
    Matrix2x2 Ainv = dagger(A);
    V              = mmul(mmul(A, V0), Ainv);
    W              = mmul(mmul(A, W0), Ainv);
}

using Word = std::vector<char>;

Matrix2x2 word_to_matrix(const Word& w) {
    Matrix2x2 U = mat_I();
    for (char g : w) {
        if (g == 'H')
            U = mmul(U, mat_H());
        else if (g == 'T')
            U = mmul(U, mat_T());
        else if (g == 't')
            U = mmul(U, mat_Tdg());
        else
            throw std::runtime_error("Unknown gate symbol in word");
    }
    return U;
}

Word invert_word(const Word& w) {
    Word inv;
    inv.reserve(w.size());
    for (auto it = w.rbegin(); it != w.rend(); ++it) {
        char g = *it;
        if (g == 'H')
            inv.push_back('H');
        else if (g == 'T')
            inv.push_back('t');
        else if (g == 't')
            inv.push_back('T');
        else
            throw std::runtime_error("Unknown gate symbol in word");
    }
    return inv;
}

std::string hash_upto_phase(const Matrix2x2& U_in) {
    Matrix2x2 U = normalize_to_su2(U_in);
    Complex   z = U[0][0];
    if (std::abs(z) < 1e-12)
        z = U[1][1];
    if (std::abs(z) > 1e-12) {
        Complex ph = z / std::abs(z);
        Complex s  = std::conj(ph);
        for (int r = 0; r < 2; r++)
            for (int c = 0; c < 2; c++)
                U[r][c] *= s;
    }
    auto        rnd = [](double x) -> long long { return (long long)std::llround(x * 1e10); };
    std::string k;
    k.reserve(128);
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) {
            k += std::to_string(rnd(std::real(U[r][c]))) + ",";
            k += std::to_string(rnd(std::imag(U[r][c]))) + ";";
        }
    return k;
}

struct NetElem {
    Matrix2x2 U;
    Word      w;
};

std::vector<NetElem> build_net(int max_len) {
    std::vector<NetElem> net;
    net.reserve(20000);

    std::unordered_map<std::string, int> seen;
    auto                                 push = [&](const Matrix2x2& U, const Word& w) {
        std::string key = hash_upto_phase(U);
        if (seen.find(key) != seen.end())
            return;
        seen.emplace(key, (int)net.size());
        net.push_back({U, w});
    };

    std::queue<NetElem> q;
    push(mat_I(), {});
    q.push({mat_I(), {}});

    const std::array<std::pair<char, Matrix2x2>, 3> gens = {{
        {'H', mat_H()},
        {'T', mat_T()},
        {'t', mat_Tdg()},
    }};

    while (!q.empty()) {
        NetElem cur = q.front();
        q.pop();
        if ((int)cur.w.size() == max_len)
            continue;
        for (auto& [sym, G] : gens) {
            NetElem nxt;
            nxt.U = mmul(cur.U, G);
            nxt.w = cur.w;
            nxt.w.push_back(sym);

            std::string key = hash_upto_phase(nxt.U);
            if (seen.find(key) != seen.end())
                continue;
            seen.emplace(key, (int)net.size());
            net.push_back(nxt);
            q.push(nxt);
        }
    }
    return net;
}

const NetElem& nearest_in_net(const std::vector<NetElem>& net, const Matrix2x2& U) {
    double best   = std::numeric_limits<double>::infinity();
    size_t best_i = 0;
    for (size_t i = 0; i < net.size(); i++) {
        double d = dist_phase_invariant(net[i].U, U);
        if (d < best) {
            best   = d;
            best_i = i;
        }
    }
    return net[best_i];
}

Word sk_synthesize(const Matrix2x2& U_in, int depth, const std::vector<NetElem>& net) {
    Matrix2x2 U = normalize_to_su2(U_in);

    if (depth <= 0) {
        return nearest_in_net(net, U).w;
    }

    const NetElem& U0e = nearest_in_net(net, U);
    Matrix2x2      U0  = U0e.U;
    Word           w0  = U0e.w;

    Matrix2x2 Delta = mmul(U, dagger(U0));
    Delta           = normalize_to_su2(Delta);

    Matrix2x2 V, W;
    decompose_near_identity_commutator(Delta, V, W);

    Word wV = sk_synthesize(V, depth - 1, net);
    Word wW = sk_synthesize(W, depth - 1, net);

    Word wVinv = invert_word(wV);
    Word wWinv = invert_word(wW);

    Word out;
    out.reserve(wV.size() + wW.size() + wVinv.size() + wWinv.size() + w0.size());
    out.insert(out.end(), wV.begin(), wV.end());
    out.insert(out.end(), wW.begin(), wW.end());
    out.insert(out.end(), wVinv.begin(), wVinv.end());
    out.insert(out.end(), wWinv.begin(), wWinv.end());
    out.insert(out.end(), w0.begin(), w0.end());
    return out;
}

} // namespace

std::vector<std::shared_ptr<QGate>> approximate_ry_sk(uint32_t target, double theta, double eps) {
    if (!(eps > 0.0) || !std::isfinite(eps))
        throw std::invalid_argument("approximate_ry_sk: eps must be finite and > 0");
    if (!std::isfinite(theta))
        throw std::invalid_argument("approximate_ry_sk: theta must be finite");

    Matrix2x2 U = mat_Ry(theta);
    U           = normalize_to_su2(U);

    static const int                  NET_LEN = 8;
    static const std::vector<NetElem> NET     = build_net(NET_LEN);

    const int MAX_DEPTH = 5;

    Word   best_word;
    double best_d = std::numeric_limits<double>::infinity();

    for (int depth = 0; depth <= MAX_DEPTH; depth++) {
        Word      w  = sk_synthesize(U, depth, NET);
        Matrix2x2 Uw = normalize_to_su2(word_to_matrix(w));
        double    d  = dist_phase_invariant(Uw, U);

        if (d < best_d) {
            best_d    = d;
            best_word = std::move(w);
        }
        if (d <= eps)
            break;
    }

    std::vector<std::shared_ptr<QGate>> gates;
    gates.reserve(best_word.size());
    for (char g : best_word) {
        if (g == 'H')
            gates.push_back(std::make_shared<H>(target));
        else if (g == 'T')
            gates.push_back(std::make_shared<T>(target));
        else if (g == 't')
            gates.push_back(std::make_shared<Tdg>(target));
        else
            throw std::runtime_error("approximate_ry_sk: unknown symbol in synthesized word");
    }
    return gates;
}

} // namespace xyz
