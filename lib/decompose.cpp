#include "qgate.hpp"
#include "transpile.hpp"

#include <cmath>
#include <vector>

namespace xyz {
namespace {

std::vector<double> find_thetas(const std::vector<double>& alphas) {
    uint32_t                         size = alphas.size();
    std::vector<double>              thetas(size, 0.0);
    std::vector<std::vector<double>> mat(size, std::vector<double>(size, 0.0));

    for (uint32_t i = 0; i < size; i++) {
        for (uint32_t j = 0; j < size; j++) {
            uint32_t gray_j = j ^ (j >> 1);
            uint32_t pop    = __builtin_popcount(i & gray_j);
            mat[i][j]       = (pop & 1) ? -1.0 : 1.0;
        }
    }

    for (uint32_t i = 0; i < size; i++) {
        uint32_t pivot = i;
        for (uint32_t j = i + 1; j < size; j++)
            if (std::abs(mat[j][i]) > std::abs(mat[pivot][i]))
                pivot = j;
        if (pivot != i) {
            std::swap(mat[i], mat[pivot]);
            std::swap(thetas[i], thetas[pivot]);
        }
        for (uint32_t j = i + 1; j < size; j++) {
            double factor = mat[j][i] / mat[i][i];
            for (uint32_t k = i; k < size; k++)
                mat[j][k] -= factor * mat[i][k];
        }
    }

    for (uint32_t i = 0; i < size; i++)
        thetas[i] = alphas[i];

    for (uint32_t i = 0; i < size; i++) {
        for (uint32_t j = i + 1; j < size; j++) {
            double factor = mat[i][j] / mat[i][i];
            for (uint32_t k = j; k < size; k++)
                mat[i][k] -= factor * mat[j][k];
            thetas[i] -= factor * thetas[j];
        }
    }

    for (uint32_t i = 0; i < size; i++)
        thetas[i] /= mat[i][i];

    return thetas;
}

} // namespace

std::vector<std::shared_ptr<QGate>> decompose_mcry(const MCRY& gate) {
    uint32_t            num_controls = gate.ctrls.size();
    uint32_t            table_size   = 1 << num_controls;
    std::vector<double> rotation_table(table_size, 0.0);

    uint32_t rotated_index = 0;
    for (uint32_t i = 0; i < num_controls; i++)
        if (gate.phases[i])
            rotated_index += 1 << i;

    rotation_table[rotated_index] = gate.theta;

    std::vector<double> thetas = find_thetas(rotation_table);

    std::vector<std::shared_ptr<QGate>> gates;
    uint32_t                            prev_gray = 0;

    for (uint32_t i = 0; i < table_size; i++) {
        uint32_t curr_gray = (i + 1) ^ ((i + 1) >> 1);
        if (i == table_size - 1)
            curr_gray = 0;
        uint32_t diff       = curr_gray ^ prev_gray;
        uint32_t control_id = __builtin_ctz(diff);
        prev_gray           = curr_gray;

        gates.push_back(std::make_shared<RY>(gate.target, thetas[i]));
        gates.push_back(std::make_shared<CX>(gate.ctrls[control_id], true, gate.target));
    }

    return gates;
}

} // namespace xyz
