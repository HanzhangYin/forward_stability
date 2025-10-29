# Script to compute and plot E[FS(u,v)] vs N
# This runs forward stability computations for different values of N
# and plots how the expected value changes with N

from sage.all import *
from collections import Counter
import numpy as np
import time
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import matplotlib.pyplot as plt

# --- Define constants ---
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)

# --- Configuration ---
N_START = 20
N_END = 90
N_STEP = 5  # Compute every 5th value (30, 35, 40, ..., 80)
NUM_RANDOM_SAMPLES = 100000  # Number of samples per N
POOL_SIZE = 7000

# =============================================================================
# HELPER FUNCTIONS (copied from forward_stability_distribution.py)
# =============================================================================

def FS_ambient(u, v):
    """
    Computes the Forward Stability Number using the ambient formula.
    This is faster than integer_support_from_formula (O(n) vs O(n²)).
    
    Args:
        u, v: Permutation objects
        
    Returns:
        int: The Forward Stability Number FS(u,v)
    """
    n = len(u)
    du = u.to_lehmer_cocode()
    dv = v.to_lehmer_cocode()
    nz_u = [1 if x > 0 else 0 for x in du]
    nz_v = [1 if x > 0 else 0 for x in dv]
    tail_u = [0] * (n + 2)
    tail_v = [0] * (n + 2)
    for i in range(n, 0, -1):
        tail_u[i] = tail_u[i + 1] + nz_u[i - 1]
        tail_v[i] = tail_v[i + 1] + nz_v[i - 1]
    raw = max(tail_u[i] + tail_v[i] + i - 1 for i in range(1, n + 1))
    return max(n, raw)

def process_pair(pair):
    """Worker function for multiprocessing."""
    try:
        u, v = pair
        return FS_ambient(u, v)
    except Exception as e:
        print(f"Error processing pair {pair}: {e}")
        return 0

def compute_expectation(frequencies):
    """Compute expected value from frequency distribution."""
    total_count = sum(frequencies.values())
    if total_count == 0:
        return 0.0
    weighted_sum = sum(value * count for value, count in frequencies.items())
    return weighted_sum / total_count

# =============================================================================
# MAIN COMPUTATION FUNCTION
# =============================================================================

def compute_expectation_for_n(n, num_samples=NUM_RANDOM_SAMPLES):
    """
    Computes E[FS(u,v)] for a given value of n.
    
    Args:
        n: The size of the symmetric group S_n
        num_samples: Number of random pairs to sample
        
    Returns:
        float: The expected value E[FS(u,v)]
    """
    print(f"\n{'='*60}")
    print(f"Computing for S_{n}")
    print(f"{'='*60}")
    
    perms_obj = Permutations(n)
    
    # Generate random pairs
    print(f"Pre-generating pool of {min(POOL_SIZE, num_samples * 2)} permutations...")
    pool_size = min(POOL_SIZE, num_samples * 2)
    perm_pool = [perms_obj.random_element() for _ in range(pool_size)]
    print(f"Pool generation complete.")
    
    print(f"Generating {num_samples} random index pairs...")
    indices_u = np.random.randint(0, pool_size, num_samples)
    indices_v = np.random.randint(0, pool_size, num_samples)
    
    print("Building task list from random pairs...")
    task_list = [(perm_pool[indices_u[i]], perm_pool[indices_v[i]]) for i in range(num_samples)]
    print(f"Task list of {len(task_list)} pairs built.")
    
    # Parallel processing
    num_cores = cpu_count()
    print(f"\nStarting parallel processing with {num_cores} cores...")
    start_time = time.time()
    
    with Pool(processes=num_cores) as pool:
        chunksize = max(1, len(task_list) // (num_cores * 4))
        results = list(tqdm(
            pool.imap(process_pair, task_list, chunksize=chunksize),
            total=len(task_list),
            desc=f"Processing S_{n}",
            unit="pair"
        ))
    
    end_time = time.time()
    print(f"\nParallel processing finished in {end_time - start_time:.4f} seconds.")
    
    # Compute expectation
    support_frequencies = Counter(results)
    expectation = compute_expectation(support_frequencies)
    
    print(f"E[FS(u,v)] for S_{n} = {expectation:.4f}")
    
    return expectation

# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    """Main function to compute and plot E[FS(u,v)] vs N."""
    
    print("="*60)
    print("EXPECTATION VALUE vs N ANALYSIS")
    print("="*60)
    print(f"Range: N = {N_START} to {N_END} (step {N_STEP})")
    print(f"Samples per N: {NUM_RANDOM_SAMPLES}")
    print("="*60)
    
    # Compute expectations for each N
    n_values = list(range(N_START, N_END + 1, N_STEP))
    expectations = []
    
    for n in n_values:
        exp_val = compute_expectation_for_n(n, NUM_RANDOM_SAMPLES)
        expectations.append(exp_val)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY RESULTS")
    print("="*60)
    for n, exp in zip(n_values, expectations):
        print(f"S_{n}: E[FS(u,v)] = {exp:.4f}")
    print("="*60)
    
    # Plot the results
    print("\nGenerating plot...")
    plt.figure(figsize=(12, 7))
    plt.plot(n_values, expectations, 'o-', linewidth=2, markersize=8, color='steelblue')
    plt.xlabel('N (Symmetric Group S_N)', fontsize=14)
    plt.ylabel('E[FS(u,v)]', fontsize=14)
    plt.title(f'Expected Forward Stability Number vs N\n({NUM_RANDOM_SAMPLES:,} samples per N)', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.xticks(n_values, rotation=45)
    
    # Add data point labels
    for n, exp in zip(n_values, expectations):
        plt.annotate(f'{exp:.2f}', 
                    xy=(n, exp), 
                    xytext=(0, 10), 
                    textcoords='offset points',
                    ha='center',
                    fontsize=9,
                    alpha=0.7)
    
    plt.tight_layout()
    filename = f'expectation_vs_n_{N_START}_to_{N_END}.png'
    plt.savefig(filename, dpi=150)
    print(f"Plot saved as '{filename}'")
    
    # Also save data to text file
    data_filename = f'expectation_data_{N_START}_to_{N_END}.txt'
    with open(data_filename, 'w') as f:
        f.write("N\tE[FS(u,v)]\n")
        for n, exp in zip(n_values, expectations):
            f.write(f"{n}\t{exp:.6f}\n")
    print(f"Data saved as '{data_filename}'")
    
    print("\nScript finished!")

if __name__ == "__main__":
    main()

