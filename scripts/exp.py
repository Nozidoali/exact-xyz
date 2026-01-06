import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import gmean

from xyz_wrapper import prepare_xyz_circuit
from qualtran_wrapper import prepare_qualtran_circuit


def generate_sparse_state(n, cardinality, seed, uniform=True):
    np.random.seed(seed)
    dim = 2**n
    cardinality = min(cardinality, dim)
    
    indices = np.random.choice(dim, size=cardinality, replace=False)
    coeffs = np.zeros(dim)
    
    if uniform:
        # Uniform coefficients: all non-zero entries have equal magnitude
        coeffs[indices] = 1.0
    else:
        # Random coefficients: sampled from normal distribution
        coeffs[indices] = np.random.randn(cardinality)
    
    coeffs = coeffs / np.linalg.norm(coeffs)
    
    return coeffs


def run_benchmark(n_values, card_values, eps_values, num_repeats=6, seed=42, uniform=True):
    results = []
    total = sum(1 for n in n_values for c in card_values if c <= 2**n) * len(eps_values) * num_repeats
    count = 0
    
    for n in n_values:
        for card in card_values:
            if card > 2**n:
                continue
            for eps in eps_values:
                for repeat in range(num_repeats):
                    count += 1
                    # Use different seed for each repeat to get different random states
                    current_seed = seed + repeat
                    coeffs = generate_sparse_state(n, card, current_seed, uniform=uniform)
                    
                    try:
                        xyz_result = prepare_xyz_circuit(coeffs, eps)
                        
                        result = {
                            'n': n,
                            'cardinality': card,
                            'eps': eps,
                            'repeat': repeat,
                            'xyz_t': xyz_result['t'],
                            'xyz_gates': xyz_result['total_gates'],
                            'xyz_qubits': xyz_result['qubits'],
                        }
                        
                        qt_result = prepare_qualtran_circuit(coeffs, eps)
                        result.update({
                            'qt_t': qt_result['t'],
                            'qt_clifford': qt_result['clifford'],
                            'qt_qubits': qt_result['qubits'],
                            'qt_phase_bits': qt_result['phase_bitsize'],
                        })
                        
                        results.append(result)
                        print(f"[{count}/{total}] n={n}, card={card}, eps={eps:.0e}, rep={repeat+1}/{num_repeats}: XYZ T={result['xyz_t']:3d}", end="")
                        print(f", QT T={result['qt_t']:4d}")
                            
                    except Exception as e:
                        print(f"[{count}/{total}] Error n={n}, card={card}, eps={eps:.0e}, rep={repeat+1}: {e}")
    
    return pd.DataFrame(results)


def main():
    print("=== XYZ vs Qualtran Benchmark ===\n")
    
    n_values = [3, 4, 5]
    card_values = [2, 4, 8, 16]
    eps_values = [1e-1, 1e-2, 1e-3]
    num_repeats = 6
    uniform = True  # Use uniform coefficients by default
    
    print(f"Config: n={n_values}, card={card_values}, eps={eps_values}")
    print(f"Uniform coefficients: {uniform}")
    print(f"Number of repeats: {num_repeats}\n")
    
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    
    # Run benchmarks for each epsilon and save separately
    for eps in eps_values:
        print(f"\n{'='*60}")
        print(f"Running benchmarks for eps={eps:.0e}")
        print(f"{'='*60}\n")
        
        df = run_benchmark(n_values, card_values, [eps], num_repeats=num_repeats, seed=42, uniform=uniform)
        
        if df.empty:
            print(f"No results for eps={eps:.0e}")
            continue
        
        # Calculate statistics for each configuration
        grouped = df.groupby(['n', 'cardinality', 'eps'])
        stats_list = []
        
        for (n, card, eps_val), group in grouped:
            stats = {
                'n': n,
                'cardinality': card,
                'eps': eps_val,
                'xyz_t_mean': group['xyz_t'].mean(),
                'xyz_t_std': group['xyz_t'].std(),
                'xyz_gates_mean': group['xyz_gates'].mean(),
                'xyz_gates_std': group['xyz_gates'].std(),
                'xyz_qubits_mean': group['xyz_qubits'].mean(),
                'xyz_qubits_std': group['xyz_qubits'].std(),
                'qt_t_mean': group['qt_t'].mean(),
                'qt_t_std': group['qt_t'].std(),
                'qt_clifford_mean': group['qt_clifford'].mean(),
                'qt_clifford_std': group['qt_clifford'].std(),
                'qt_qubits_mean': group['qt_qubits'].mean(),
                'qt_qubits_std': group['qt_qubits'].std(),
                'qt_phase_bits_mean': group['qt_phase_bits'].mean(),
                'qt_phase_bits_std': group['qt_phase_bits'].std(),
            }
            stats_list.append(stats)
        
        stats_df = pd.DataFrame(stats_list)
        
        # Calculate geometric mean of T-count ratios
        stats_df['t_ratio_mean'] = stats_df['xyz_t_mean'] / stats_df['qt_t_mean'].replace(0, np.nan)
        geomean_ratio = gmean(stats_df['t_ratio_mean'].dropna())
        
        # Save raw results
        raw_csv_file = output_dir / f"benchmark_results_eps{eps:.0e}_raw.csv"
        df.to_csv(raw_csv_file, index=False)
        print(f"\nSaved raw results to {raw_csv_file}")
        
        # Save aggregated statistics
        stats_csv_file = output_dir / f"benchmark_results_eps{eps:.0e}.csv"
        stats_df.to_csv(stats_csv_file, index=False)
        print(f"Saved statistics to {stats_csv_file}")
        
        print(f"\n=== Summary for eps={eps:.0e} ===")
        print(f"XYZ avg T: {stats_df['xyz_t_mean'].mean():.1f} ± {stats_df['xyz_t_std'].mean():.1f}")
        print(f"QT avg T:  {stats_df['qt_t_mean'].mean():.1f} ± {stats_df['qt_t_std'].mean():.1f}")
        print(f"Geometric mean T-count ratio (XYZ/QT): {geomean_ratio:.3f}x")
    
    print("\n" + "="*60)
    print("✓ Complete! Run 'python plot.py' to generate plots.")
    print("="*60)


if __name__ == "__main__":
    main()
