#include "qgate.hpp"
#include "transpile.hpp"

#include <cmath>
#include <vector>

namespace xyz {
namespace {

std::vector<double> find_thetas(const std::vector<double>& alphas) {
    uint32_t             size = alphas.size();
    std::vector<double>  thetas(size, 0.0);
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
    uint32_t             num_controls = gate.ctrls.size();
    uint32_t             table_size   = 1 << num_controls;
    std::vector<double>  rotation_table(table_size, 0.0);
    
    uint32_t rotated_index = 0;
    for (uint32_t i = 0; i < num_controls; i++)
        if (gate.phases[i])
            rotated_index += 1 << i;
    
    rotation_table[rotated_index] = gate.theta;
    
    std::vector<double> thetas = find_thetas(rotation_table);
    
    std::vector<std::shared_ptr<QGate>> gates;
    uint32_t                            prev_gray = 0;
    
    for (uint32_t i = 0; i < table_size; i++) {
        uint32_t curr_gray  = (i + 1) ^ ((i + 1) >> 1);
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

std::vector<std::shared_ptr<QGate>> decompose_mcry_qrom(const MCRY& gate, double eps) {
    uint32_t             num_controls = gate.ctrls.size();
    uint32_t             table_size   = 1 << num_controls;
    std::vector<double>  rotation_table(table_size, 0.0);
    
    uint32_t rotated_index = 0;
    for (uint32_t i = 0; i < num_controls; i++)
        if (gate.phases[i])
            rotated_index += 1 << i;
    
    rotation_table[rotated_index] = gate.theta;
    
    std::vector<std::shared_ptr<QGate>> gates;
    uint32_t                            num_ancilla = num_controls;
    uint32_t                            total_qubits = gate.target + 1;
    for (uint32_t c : gate.ctrls)
        total_qubits = std::max(total_qubits, c + 1);
    uint32_t ancilla_start = total_qubits;
    
    for (uint32_t i = 0; i < num_ancilla; i++)
        gates.push_back(std::make_shared<H>(ancilla_start + i));
    
    for (uint32_t addr = 0; addr < table_size; addr++) {
        if (std::abs(rotation_table[addr]) < 1e-10)
            continue;
        
        std::vector<uint32_t> qrom_ctrls;
        for (uint32_t i = 0; i < num_controls; i++)
            qrom_ctrls.push_back(ancilla_start + i);
        qrom_ctrls.insert(qrom_ctrls.end(), gate.ctrls.begin(), gate.ctrls.end());
        
        std::vector<bool> qrom_phases;
        for (uint32_t i = 0; i < num_controls; i++)
            qrom_phases.push_back((addr >> i) & 1);
        qrom_phases.insert(qrom_phases.end(), gate.phases.begin(), gate.phases.end());
        
        auto rz_word = sk::synthesize_rz(rotation_table[addr], eps);
        for (auto g : rz_word) {
            if (g == sk::Gate::H) {
                gates.push_back(std::make_shared<H>(gate.target));
            } else if (g == sk::Gate::T) {
                std::vector<uint32_t> t_ctrls = qrom_ctrls;
                std::vector<bool>     t_phases = qrom_phases;
                gates.push_back(std::make_shared<MCRY>(t_ctrls, t_phases, M_PI / 4.0, gate.target));
            } else {
                std::vector<uint32_t> t_ctrls = qrom_ctrls;
                std::vector<bool>     t_phases = qrom_phases;
                gates.push_back(std::make_shared<MCRY>(t_ctrls, t_phases, -M_PI / 4.0, gate.target));
            }
        }
    }
    
    for (int i = num_ancilla - 1; i >= 0; i--)
        gates.push_back(std::make_shared<H>(ancilla_start + i));
    
    return gates;
}

} // namespace xyz

