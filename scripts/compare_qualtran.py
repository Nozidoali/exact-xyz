#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np

import attrs
if not hasattr(attrs, '__version_info__') or attrs.__version__ >= '22.2':
    import sys
    print(f"Warning: attrs version {attrs.__version__} may have compatibility issues")
    
from qualtran.bloqs.state_preparation.state_preparation_via_rotation import StatePreparationViaRotations
from qualtran.bloqs.rotations.phase_gradient import PhaseGradientState
from qualtran import BloqBuilder
from qualtran.resource_counting import get_cost_value, QECGatesCost


def generate_random_state(state_bitsize, seed=None):
    if seed is not None:
        np.random.seed(seed)
    state_coefs = np.random.randn(2**state_bitsize) + 1j*np.random.randn(2**state_bitsize)
    return state_coefs / np.linalg.norm(state_coefs)


def build_qualtran_circuit(state_coefs, phase_bitsize=8):
    state_bitsize = int(np.log2(len(state_coefs)))
    qsp = StatePreparationViaRotations(phase_bitsize=phase_bitsize, state_coefficients=tuple(state_coefs))
    return qsp


def get_qualtran_costs(qsp):
    try:
        t_complexity = qsp.t_complexity()
        return {
            't': t_complexity.t,
            'clifford': t_complexity.clifford,
            'rotation': getattr(t_complexity, 'rotations', 0),
        }
    except Exception as e:
        print(f"Error getting costs: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Get Qualtran state preparation QEC costs")
    parser.add_argument("-n", "--num-qubits", type=int, default=3)
    parser.add_argument("-p", "--phase-bitsize", type=int, default=8)
    parser.add_argument("-s", "--seed", type=int, default=None)
    args = parser.parse_args()
    
    print(f"State qubits: {args.num_qubits}, Phase bitsize: {args.phase_bitsize}, Seed: {args.seed}\n")
    
    state_coefs = generate_random_state(args.num_qubits, args.seed)
    qsp = build_qualtran_circuit(state_coefs, args.phase_bitsize)
    costs = get_qualtran_costs(qsp)
    
    if costs is None:
        print("Failed to get Qualtran costs")
        return
    
    print("Qualtran QEC Costs:")
    print(f"  T: {costs['t']}")
    print(f"  Clifford: {costs['clifford']}")
    print(f"  Rotation: {costs['rotation']}")


if __name__ == "__main__":
    main()
