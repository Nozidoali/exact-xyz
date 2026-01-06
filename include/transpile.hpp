#pragma once

#include "qcircuit.hpp"

#include <complex>
#include <cstdint>
#include <vector>

namespace xyz::sk {

enum class Gate : uint8_t { H, T, Tdg };

struct SU2 {
    std::complex<double> a{1, 0}, b{0, 0}, c{0, 0}, d{1, 0};
};

SU2               identity();
SU2               dagger(const SU2& u);
SU2               mul(const SU2& x, const SU2& y);
SU2               normalize_det(const SU2& u);
double            dist(const SU2& u, const SU2& v);
SU2               gate_matrix(Gate g);
SU2               rz_matrix(double theta);
SU2               word_matrix(const std::vector<Gate>& w);
std::vector<Gate> invert_word(const std::vector<Gate>& w);
std::vector<Gate> synthesize_rz(double theta, double eps);

} // namespace xyz::sk

namespace xyz {

QCircuit transpile_clifford_t(const QCircuit& in, double eps);

} // namespace xyz

