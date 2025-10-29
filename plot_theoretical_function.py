# Script to plot the theoretical function f(N) = 2N - 4/3 log_2(N)
# and optionally compare with empirical expectation values

import numpy as np
import matplotlib.pyplot as plt
import os

# --- Configuration ---
N_START = 20
N_END = 90
N_STEP = 5

def theoretical_function(n):
    """
    Compute f(n) = 2n - 4/3 log_2(n)
    
    Args:
        n: Input value
        
    Returns:
        float: Value of 2n - 4/3 log_2(n)
    """
    return 2 * n - 1.333 * np.log(n) / np.log(2)

def load_empirical_data(filename):
    """
    Load empirical expectation data from file if it exists.
    
    Args:
        filename: Path to the data file
        
    Returns:
        tuple: (n_values, expectations) or (None, None) if file doesn't exist
    """
    if not os.path.exists(filename):
        return None, None
    
    n_values = []
    expectations = []
    
    with open(filename, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                n_values.append(int(parts[0]))
                expectations.append(float(parts[1]))
    
    return n_values, expectations

def main():
    """Main function to plot theoretical function and compare with empirical data."""
    
    print("="*60)
    print("THEORETICAL FUNCTION PLOT")
    print(f"f(N) = 2N - 4/3 log_2(N)")
    print("="*60)
    print(f"Range: N = {N_START} to {N_END} (step {N_STEP})")
    print("="*60)
    
    # Generate N values
    n_values = np.arange(N_START, N_END + 1, N_STEP)
    
    # Compute theoretical function values
    theoretical_values = [theoretical_function(n) for n in n_values]
    
    # Print values
    print("\nTheoretical Values:")
    print("-" * 40)
    print(f"{'N':>5} {'2N - 4/3 log_2(N)':>15}")
    print("-" * 40)
    for n, val in zip(n_values, theoretical_values):
        print(f"{n:>5} {val:>15.4f}")
    print("-" * 40)
    
    # Try to load empirical data
    empirical_file = f'expectation_data_{N_START}_to_{N_END}.txt'
    empirical_n, empirical_exp = load_empirical_data(empirical_file)
    
    # Create plot
    print("\nGenerating plot...")
    plt.figure(figsize=(12, 7))
    
    # Plot theoretical function
    plt.plot(n_values, theoretical_values, 'o-', 
             linewidth=2, markersize=8, color='red', 
             label='f(N) = 2N - 4/3 log_2(N) (Theoretical)')
    
    # Add theoretical data point labels
    for n, val in zip(n_values, theoretical_values):
        plt.annotate(f'{val:.2f}', 
                    xy=(n, val), 
                    xytext=(0, -15), 
                    textcoords='offset points',
                    ha='center',
                    fontsize=9,
                    alpha=0.7,
                    color='red')
    
    # Plot empirical data if available
    if empirical_n is not None and empirical_exp is not None:
        print("Found empirical data - adding to plot for comparison...")
        plt.plot(empirical_n, empirical_exp, 's-', 
                linewidth=2, markersize=8, color='steelblue', 
                label='E[FS(u,v)] (Empirical)')
        
        # Add empirical data point labels
        for n, exp in zip(empirical_n, empirical_exp):
            plt.annotate(f'{exp:.2f}', 
                        xy=(n, exp), 
                        xytext=(0, 10), 
                        textcoords='offset points',
                        ha='center',
                        fontsize=9,
                        alpha=0.7,
                        color='steelblue')
        
        # Compute and print differences
        print("\nComparison (Theoretical vs Empirical):")
        print("-" * 60)
        print(f"{'N':>5} {'Theoretical':>15} {'Empirical':>15} {'Difference':>15}")
        print("-" * 60)
        for n, theo, emp in zip(empirical_n, 
                                 [theoretical_function(n) for n in empirical_n], 
                                 empirical_exp):
            diff = theo - emp
            print(f"{n:>5} {theo:>15.4f} {emp:>15.4f} {diff:>15.4f}")
        print("-" * 60)
    else:
        print(f"No empirical data found at '{empirical_file}'")
        print("Plotting theoretical function only.")
    
    # Format plot
    plt.xlabel('N (Symmetric Group S_N)', fontsize=14)
    plt.ylabel('Value', fontsize=14)
    plt.title('Theoretical Function: f(N) = 2N - 4/3 log_2(N)', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.xticks(n_values, rotation=45)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    
    # Save plot
    filename = f'theoretical_function_{N_START}_to_{N_END}.png'
    plt.savefig(filename, dpi=150)
    print(f"\nPlot saved as '{filename}'")
    
    # Save theoretical data to file
    data_filename = f'theoretical_data_{N_START}_to_{N_END}.txt'
    with open(data_filename, 'w') as f:
        f.write("N\t2N - 4/3 log_2(N)\n")
        for n, val in zip(n_values, theoretical_values):
            f.write(f"{n}\t{val:.6f}\n")
    print(f"Data saved as '{data_filename}'")
    
    print("\nScript finished!")

if __name__ == "__main__":
    main()

