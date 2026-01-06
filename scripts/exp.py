import numpy as np
import pandas as pd
from pathlib import Path

from xyz_wrapper import prepare_xyz_circuit
from qualtran_wrapper import prepare_qualtran_circuit, QUALTRAN_AVAILABLE


def generate_sparse_state(n, cardinality, seed):
    np.random.seed(seed)
    dim = 2**n
    cardinality = min(cardinality, dim)
    
    indices = np.random.choice(dim, size=cardinality, replace=False)
    coeffs = np.zeros(dim)
    coeffs[indices] = np.random.randn(cardinality)
    coeffs = coeffs / np.linalg.norm(coeffs)
    
    return coeffs


def run_benchmark(n_values, card_values, eps_values, seed=42):
    results = []
    total = sum(1 for n in n_values for c in card_values if c <= 2**n) * len(eps_values)
    count = 0
    
    for n in n_values:
        for card in card_values:
            if card > 2**n:
                continue
            for eps in eps_values:
                count += 1
                coeffs = generate_sparse_state(n, card, seed)
                
                try:
                    xyz_result = prepare_xyz_circuit(coeffs, eps)
                    
                    result = {
                        'n': n,
                        'cardinality': card,
                        'eps': eps,
                        'xyz_t': xyz_result['t'],
                        'xyz_gates': xyz_result['total_gates'],
                        'xyz_qubits': xyz_result['qubits'],
                    }
                    
                    if QUALTRAN_AVAILABLE:
                        qt_result = prepare_qualtran_circuit(coeffs, eps)
                        result.update({
                            'qt_t': qt_result['t'],
                            'qt_clifford': qt_result['clifford'],
                            'qt_qubits': qt_result['qubits'],
                            'qt_phase_bits': qt_result['phase_bitsize'],
                        })
                    
                    results.append(result)
                    print(f"[{count}/{total}] n={n}, card={card}, eps={eps:.0e}: XYZ T={result['xyz_t']:3d}", end="")
                    if QUALTRAN_AVAILABLE:
                        print(f", QT T={result['qt_t']:4d}")
                    else:
                        print()
                        
                except Exception as e:
                    print(f"[{count}/{total}] Error n={n}, card={card}: {e}")
    
    return pd.DataFrame(results)


def main():
    print("=== XYZ vs Qualtran Benchmark ===\n")
    
    n_values = [3, 4, 5, 6]
    card_values = [2, 4, 8, 16]
    eps_values = [1e-3]
    
    print(f"Config: n={n_values}, card={card_values}, eps={eps_values}\n")
    
    df = run_benchmark(n_values, card_values, eps_values, seed=42)
    
    if df.empty:
        print("No results")
        return
    
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    
    csv_file = output_dir / "benchmark_results.csv"
    df.to_csv(csv_file, index=False)
    print(f"\nSaved to {csv_file}")
    
    if QUALTRAN_AVAILABLE:
        print(f"\n=== Summary ===")
        print(f"XYZ avg T: {df['xyz_t'].mean():.1f}")
        print(f"QT avg T:  {df['qt_t'].mean():.1f}")
        df['t_ratio'] = df['xyz_t'] / df['qt_t'].replace(0, np.nan)
        print(f"Avg ratio: {df['t_ratio'].mean():.3f}x")
    
    print("\n✓ Complete! Run 'python plot.py' to generate plots.")


if __name__ == "__main__":
    main()
