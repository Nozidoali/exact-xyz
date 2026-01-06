#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

namespace xyz {

class QGate;

using Complex   = std::complex<double>;
using Matrix2x2 = std::array<std::array<Complex, 2>, 2>;

Matrix2x2 mmul(const Matrix2x2& a, const Matrix2x2& b);
Matrix2x2 dagger(const Matrix2x2& a);
Complex   tr(const Matrix2x2& a);
Complex   det2(const Matrix2x2& a);
Matrix2x2 normalize_to_su2(const Matrix2x2& u);
double    dist_phase_invariant(const Matrix2x2& U, const Matrix2x2& V);

Matrix2x2 mat_I();
Matrix2x2 mat_H();
Matrix2x2 mat_T();
Matrix2x2 mat_Tdg();
Matrix2x2 mat_Ry(double theta);

std::vector<std::shared_ptr<QGate>> approximate_ry_sk(uint32_t target, double theta, double eps);

} // namespace xyz
