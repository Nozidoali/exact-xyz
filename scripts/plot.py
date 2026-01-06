import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def plot_results(csv_file, output_dir="results"):
    df = pd.read_csv(csv_file)
    
    if df.empty:
        print("No data to plot")
        return
    
    Path(output_dir).mkdir(exist_ok=True)
    
    has_qualtran = 'qt_t' in df.columns
    
    if not has_qualtran:
        print("Qualtran data not available, plotting XYZ only")
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        x = np.arange(len(df))
        bars = ax.bar(x, df['xyz_t'])
        ax.set_xlabel('Configuration')
        ax.set_ylabel('T-count')
        ax.set_title('XYZ T-count (ours)')
        ax.set_xticks(x)
        ax.set_xticklabels([f"n={r['n']}\nc={r['cardinality']}" for _, r in df.iterrows()])
        for i, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom')
        ax.grid(True, alpha=0.3)
        
    else:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        
        # T-count comparison
        ax = axes[0]
        x = np.arange(len(df))
        width = 0.35
        bars1 = ax.bar(x - width/2, df['qt_t'], width, label='Qualtran', color='#A23B72')
        bars2 = ax.bar(x + width/2, df['xyz_t'], width, label='XYZ (ours)', color='#2E86AB')
        ax.set_xlabel('Configuration')
        ax.set_ylabel('T-count')
        ax.set_title('T-count Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels([f"n={r['n']}\nc={r['cardinality']}" for _, r in df.iterrows()])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        for bar in bars1:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=9)
        for bar in bars2:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=9)
        
        # Qubit comparison
        ax = axes[1]
        bars1 = ax.bar(x - width/2, df['qt_qubits'], width, label='Qualtran', color='#A23B72')
        bars2 = ax.bar(x + width/2, df['xyz_qubits'], width, label='XYZ (ours)', color='#2E86AB')
        ax.set_xlabel('Configuration')
        ax.set_ylabel('Total qubits')
        ax.set_title('Qubit Count Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels([f"n={r['n']}\nc={r['cardinality']}" for _, r in df.iterrows()])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        for bar in bars1:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom')
        for bar in bars2:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom')
        
        # Total gates comparison
        ax = axes[2]
        bars1 = ax.bar(x - width/2, df['qt_t'] + df['qt_clifford'], width, label='Qualtran (T+Clifford)', color='#A23B72')
        bars2 = ax.bar(x + width/2, df['xyz_gates'], width, label='XYZ (ours)', color='#2E86AB')
        ax.set_xlabel('Configuration')
        ax.set_ylabel('Total gates')
        ax.set_title('Total Gate Count Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels([f"n={r['n']}\nc={r['cardinality']}" for _, r in df.iterrows()])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        for bar in bars1:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=9)
        for bar in bars2:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    pdf_file = f'{output_dir}/comparison.pdf'
    plt.savefig(pdf_file, bbox_inches='tight')
    print(f"Saved plot to {pdf_file}")


def main():
    csv_file = Path(__file__).parent / "results" / "benchmark_results.csv"
    
    if not csv_file.exists():
        print(f"Error: {csv_file} not found. Run 'python exp.py' first.")
        return
    
    print(f"Loading data from {csv_file}")
    output_dir = Path(__file__).parent / "results"
    plot_results(csv_file, output_dir)
    print("✓ Plotting complete!")


if __name__ == "__main__":
    main()
