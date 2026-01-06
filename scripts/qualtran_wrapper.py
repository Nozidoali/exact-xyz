import numpy as np


from qualtran.bloqs.state_preparation.state_preparation_via_rotation import StatePreparationViaRotations


def prepare_qualtran_circuit(coeffs, eps=1e-3):

    phase_bitsize = max(8, int(np.ceil(-np.log2(eps / (2 * np.pi)))))
    
    qsp = StatePreparationViaRotations(
        phase_bitsize=phase_bitsize,
        state_coefficients=tuple(coeffs)
    )
    
    t_complexity = qsp.t_complexity()
    n_qubits = int(np.log2(len(coeffs)))
    
    return {
        't': t_complexity.t,
        'clifford': t_complexity.clifford,
        'rotation': getattr(t_complexity, 'rotations', 0),
        'qubits': n_qubits + phase_bitsize,
        'phase_bitsize': phase_bitsize
    }


if __name__ == "__main__":

    coeffs = np.array([0.6, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0])
    result = prepare_qualtran_circuit(coeffs)
    print(f"Qualtran: {result['qubits']} qubits, T={result['t']}, phase_bits={result['phase_bitsize']}")
