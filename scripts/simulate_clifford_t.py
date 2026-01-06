#!/usr/bin/env python3
import argparse
import numpy as np
from pathlib import Path
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


def simulate_circuit(qasm_file, verbose):
    circuit = QuantumCircuit.from_qasm_file(qasm_file)
    if verbose:
        print(f"Gates: {circuit.count_ops()}")
        print(f"{circuit}\n")
    statevector = Statevector.from_label('0' * circuit.num_qubits)
    return circuit.num_qubits, statevector.evolve(circuit).data


def print_statevector(num_qubits, state, threshold):
    print("State vector:")
    for i, amp in enumerate(state):
        if abs(amp) > threshold:
            basis = format(i, f'0{num_qubits}b')
            real, imag = amp.real, amp.imag
            if abs(imag) < threshold:
                print(f"  |{basis}⟩: {real:.6f}")
            elif abs(real) < threshold:
                print(f"  |{basis}⟩: {imag:.6f}i")
            else:
                sign = '+' if imag >= 0 else '-'
                print(f"  |{basis}⟩: {real:.6f} {sign} {abs(imag):.6f}i")


def print_probabilities(num_qubits, state, threshold):
    print("\nProbability distribution:")
    probs = np.abs(state)**2
    for i, prob in enumerate(probs):
        if prob > threshold:
            print(f"  |{format(i, f'0{num_qubits}b')}⟩: {prob:.6f}")


def main():
    parser = argparse.ArgumentParser(description="State vector simulation")
    parser.add_argument('input', type=str, help='Input QASM file or directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    parser.add_argument('-t', '--threshold', type=float, default=1e-6, help='Threshold')
    parser.add_argument('--probs', action='store_true', help='Show probabilities')
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: '{args.input}' not found")
        return 1
    
    qasm_files = list(input_path.glob('*.qasm')) if input_path.is_dir() else [input_path]
    if not qasm_files:
        print(f"Error: No .qasm files found")
        return 1
    
    for qasm_file in qasm_files:
        if len(qasm_files) > 1:
            print(f"\n{'='*60}\nProcessing: {qasm_file.name}")
        
        num_qubits, state = simulate_circuit(qasm_file, args.verbose)
        
        print(f"Simulating {num_qubits} qubit(s)\n{'='*60}")
        print_statevector(num_qubits, state, args.threshold)
        
        if args.probs:
            print_probabilities(num_qubits, state, args.threshold)
        
        print(f"\nTotal probability: {np.sum(np.abs(state)**2):.10f}")
    
    return 0


if __name__ == "__main__":
    exit(main())
