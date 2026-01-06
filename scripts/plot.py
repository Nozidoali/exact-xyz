import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import gmean


def plot_results_with_stats(csv_file, output_file):
    """Plot results with mean and standard deviation error bars."""
    df = pd.read_csv(csv_file)
    
    if df.empty:
        print(f"No data to plot in {csv_file}")
        return
    
    has_qualtran = 'qt_t_mean' in df.columns
    
    if not has_qualtran:
        print(f"Qualtran data not available in {csv_file}, skipping")
        return
    
    # Calculate geometric mean ratio
    df['t_ratio_mean'] = df['xyz_t_mean'] / df['qt_t_mean'].replace(0, np.nan)
    geomean_ratio = gmean(df['t_ratio_mean'].dropna())
    
    eps_val = df['eps'].iloc[0]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    x = np.arange(len(df))
    width = 0.35
    
    # T-count comparison with error bars
    ax = axes[0]
    bars1 = ax.bar(x - width/2, df['qt_t_mean'], width, 
                   yerr=df['qt_t_std'], label='Qualtran', 
                   color='#A23B72', capsize=5, alpha=0.8)
    bars2 = ax.bar(x + width/2, df['xyz_t_mean'], width, 
                   yerr=df['xyz_t_std'], label='XYZ (ours)', 
                   color='#2E86AB', capsize=5, alpha=0.8)
    ax.set_xlabel('Configuration', fontsize=12)
    ax.set_ylabel('T-count', fontsize=12)
    ax.set_title(f'T-count Comparison (ε={eps_val:.0e})', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"n={int(r['n'])}\nc={int(r['cardinality'])}" for _, r in df.iterrows()], fontsize=10)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height)}', ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Qubit comparison with error bars
    ax = axes[1]
    bars1 = ax.bar(x - width/2, df['qt_qubits_mean'], width, 
                   yerr=df['qt_qubits_std'], label='Qualtran', 
                   color='#A23B72', capsize=5, alpha=0.8)
    bars2 = ax.bar(x + width/2, df['xyz_qubits_mean'], width, 
                   yerr=df['xyz_qubits_std'], label='XYZ (ours)', 
                   color='#2E86AB', capsize=5, alpha=0.8)
    ax.set_xlabel('Configuration', fontsize=12)
    ax.set_ylabel('Total qubits', fontsize=12)
    ax.set_title(f'Qubit Count Comparison (ε={eps_val:.0e})', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"n={int(r['n'])}\nc={int(r['cardinality'])}" for _, r in df.iterrows()], fontsize=10)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    
    for bar in bars1:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height)}', ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Total gates comparison with error bars
    ax = axes[2]
    qt_total_mean = df['qt_t_mean'] + df['qt_clifford_mean']
    qt_total_std = np.sqrt(df['qt_t_std']**2 + df['qt_clifford_std']**2)
    
    bars1 = ax.bar(x - width/2, qt_total_mean, width, 
                   yerr=qt_total_std, label='Qualtran (T+Clifford)', 
                   color='#A23B72', capsize=5, alpha=0.8)
    bars2 = ax.bar(x + width/2, df['xyz_gates_mean'], width, 
                   yerr=df['xyz_gates_std'], label='XYZ (ours)', 
                   color='#2E86AB', capsize=5, alpha=0.8)
    ax.set_xlabel('Configuration', fontsize=12)
    ax.set_ylabel('Total gates', fontsize=12)
    ax.set_title(f'Total Gate Count (ε={eps_val:.0e})', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"n={int(r['n'])}\nc={int(r['cardinality'])}" for _, r in df.iterrows()], fontsize=10)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    
    for bar in bars1:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height)}', ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Add geometric mean text
    fig.text(0.5, 0.02, f'Geometric Mean T-count Ratio (XYZ/Qualtran): {geomean_ratio:.3f}x', 
             ha='center', fontsize=13, fontweight='bold', 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"Saved plot to {output_file}")
    plt.close()


def main():
    results_dir = Path(__file__).parent / "results"
    
    if not results_dir.exists():
        print(f"Error: {results_dir} not found. Run 'python exp.py' first.")
        return
    
    # Look for statistics files for each epsilon
    eps_values = [1e-1, 1e-2, 1e-3]
    
    print("=== Generating plots for all epsilon values ===\n")
    
    found_any = False
    for eps in eps_values:
        csv_file = results_dir / f"benchmark_results_eps{eps:.0e}.csv"
        
        if not csv_file.exists():
            print(f"Warning: {csv_file} not found, skipping...")
            continue
        
        found_any = True
        print(f"Processing eps={eps:.0e}...")
        
        output_file = results_dir / f"comparison_eps{eps:.0e}.pdf"
        plot_results_with_stats(csv_file, output_file)
    
    if not found_any:
        print("\nNo benchmark result files found. Run 'python exp.py' first.")
        return
    
    print("\n✓ Plotting complete!")


if __name__ == "__main__":
    main()
