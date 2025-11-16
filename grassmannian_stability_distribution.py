# Import necessary libraries
from sage.all import *
from collections import Counter
import random
import numpy as np
import time  # For timing the parallel processing
from multiprocessing import Pool, cpu_count  # <-- The new imports
from tqdm import tqdm  # For progress bar
from itertools import combinations
import pickle
import os

# --- Define constants ---
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)

# --- Configuration ---
# Change this value to compute for different Grassmannian permutations in S_n
N = 30

# Random sampling configuration
# Set to None to compute ALL pairs, or set to a number to randomly sample that many pairs
NUM_RANDOM_SAMPLES = 50000  # e.g., NUM_RANDOM_SAMPLES random pairs instead of all pairs

# Pool size for pre-generated permutations (for faster sampling)
POOL_SIZE = 10000  # Medium pool size for balance of speed and diversity

# Cache directory for pre-computed Grassmannian permutations
CACHE_DIR = ".grassmannian_cache"

# =============================================================================
# GRASSMANNIAN PERMUTATION GENERATION (OPTIMIZED)
# =============================================================================

def grassmannian_perms_direct(n):
    """
    FAST direct generation of Grassmannian permutations using combinatorial structure.

    A Grassmannian permutation has at most one descent. We generate them by:
    1. Identity permutation (no descents)
    2. For each descent position k (1 to n-1):
       - Choose k elements for positions 1..k (in increasing order)
       - Remaining elements go in positions k+1..n (in increasing order)
       - Require: max(first k elements) > min(last n-k elements) for descent

    This generates exactly 2^n - n permutations.

    Args:
        n: Size of the symmetric group

    Returns:
        list: List of Grassmannian permutations in S_n
    """
    P = Permutations(n)
    result = []
    expected_count = 2**n - n

    # Identity permutation (no descent)
    result.append(P(list(range(1, n+1))))

    # Setup progress bar with better visibility
    pbar = tqdm(
        total=expected_count,
        desc="  Generating Grassmannian perms",
        unit=" perm",
        unit_scale=False,
        bar_format='{desc}: {n_fmt}/{total_fmt} [{percentage:3.0f}%] {rate_fmt} {elapsed}<{remaining}',
        mininterval=0.1,  # Update display at least every 0.1 seconds
        miniters=1  # Update on every iteration for small n
    )
    pbar.update(1)  # Count the identity permutation

    # For each descent position k
    for k in range(1, n):
        # Choose which k elements appear in first k positions
        for first_subset in combinations(range(1, n+1), k):
            first_part = list(first_subset)  # Already sorted from combinations
            second_part = sorted(set(range(1, n+1)) - set(first_part))

            # Check if this creates a descent: max(first) > min(second)
            if first_part[-1] > second_part[0]:
                perm_list = first_part + second_part
                result.append(P(perm_list))
                pbar.update(1)

    pbar.close()
    print()  # Add newline after progress bar for cleaner output
    return result

def load_cached_grassmannian(n):
    """
    Load cached Grassmannian permutations from disk if available.

    Args:
        n: Size of the symmetric group

    Returns:
        list or None: Cached permutations if available, None otherwise
    """
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    cache_file = os.path.join(CACHE_DIR, f"grass_n{n}.pkl")

    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            print(f"Warning: Could not load cache: {e}")
            return None
    return None

def save_cached_grassmannian(n, perms):
    """
    Save generated Grassmannian permutations to disk for future use.

    Args:
        n: Size of the symmetric group
        perms: List of permutations to cache
    """
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    cache_file = os.path.join(CACHE_DIR, f"grass_n{n}.pkl")

    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(perms, f)
    except Exception as e:
        print(f"Warning: Could not save cache: {e}")

def get_grassmannian_perms(n, use_cache=True):
    """
    Get Grassmannian permutations with caching support.

    Args:
        n: Size of the symmetric group
        use_cache: Whether to use caching (default True)

    Returns:
        list: List of Grassmannian permutations in S_n
    """
    expected = 2**n - n

    # Try to load from cache first
    if use_cache:
        print(f"  Checking cache for n={n}...")
        cached = load_cached_grassmannian(n)
        if cached is not None:
            print(f"  ✓ Loaded {len(cached)} Grassmannian permutations from cache (instant!)")
            print(f"     Note: No progress bar shown when loading from cache")
            return cached
        else:
            print(f"  Cache miss - will generate from scratch")

    # Generate directly (fast!)
    print(f"  Generating {expected:,} Grassmannian permutations using direct construction...")
    print(f"     Progress bar will appear below:")
    start = time.time()
    perms = grassmannian_perms_direct(n)
    elapsed = time.time() - start
    print(f"  ✓ Generated {len(perms):,} permutations in {elapsed:.3f} seconds ({len(perms)/elapsed:.1f} perms/sec)")

    # Save to cache for future use
    if use_cache:
        save_cached_grassmannian(n, perms)
        print(f"  ✓ Saved to cache for future use")

    return perms

# =============================================================================
# TOP-LEVEL WORKER FUNCTION (FOR MULTIPROCESSING)
# =============================================================================
# This function MUST be at the top level so it can be "pickled"
# and sent to child processes. It depends on the helper functions below.
#
def process_pair(pair):
    """
    Worker function that processes a single pair (u, v).
    This is what each core will run in parallel.

    Args:
        pair (tuple): A tuple (u, v) of two Permutation objects.

    Returns:
        int: The Forward Stability Number FS(u,v).
    """
    try:
        u, v = pair
        return FS_ambient(u, v)  # Returns the value directly
    except Exception as e:
        print(f"Error processing pair {pair}: {e}")
        return 0 # Return a default value on error


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
# These must be defined before the worker function that uses them.

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

# =============================================================================
# OLD FUNCTIONS (kept for reference, not used by default)
# =============================================================================

def integer_support(w):
    """
    Finds the set of indices i for which w is not in S_i.
    This corresponds to {i | FS(w) > i}.
    """
    to_return = set()
    n = len(w)
    for i in range(_sage_const_1, n):
        for j in range(_sage_const_1, i + _sage_const_1):
            if w[j - _sage_const_1] > i:
                to_return.add(i)
                break # Found one, no need to check further for this i
    return to_return

def precompute_lambda_lookups(w, max_k):
    """
    Precomputes the cumulative sums of the theta (Lehmer)
    and dual_theta (Dual Lehmer) arrays for a permutation w.
    """
    n = len(w)
    lehmer_code = w.to_lehmer_code()
    dual_lehmer_code = w.to_lehmer_cocode()
    cum_theta = [0] * (max_k + 1)
    cum_dual_theta = [0] * (max_k + 1)

    for k in range(1, n + 1):
        cum_theta[k] = cum_theta[k - 1]
        cum_dual_theta[k] = cum_dual_theta[k - 1]
        if lehmer_code[k - 1] > 0:
            cum_theta[k] += 1
        if dual_lehmer_code[k - 1] > 0:
            cum_dual_theta[k] += 1

    for k in range(n + 1, max_k + 1):
        cum_theta[k] = cum_theta[n]
        cum_dual_theta[k] = cum_dual_theta[n]

    return cum_theta, cum_dual_theta

def integer_support_from_formula(u, v):
    """
    Calculates the set of unstable positions for the product S_u * S_v
    using fast O(1) lookups based on precomputation.
    """
    to_return = set()
    n = len(u)
    max_n = 2 * n + _sage_const_1

    u_cum_theta, u_cum_dual = precompute_lambda_lookups(u, max_n)
    v_cum_theta, v_cum_dual = precompute_lambda_lookups(v, max_n)

    to_return.update(integer_support(u))
    to_return.update(integer_support(v))

    for j in range(_sage_const_1, max_n):
        if j in to_return:
            continue

        for i in range(_sage_const_1, max_n):
           start, end = min(i, j), max(i, j)
           ord_lambda_u = u_cum_theta[end] - u_cum_theta[start - 1]
           ord_lambda_v = v_cum_theta[end] - v_cum_theta[start - 1]

           if ord_lambda_u + ord_lambda_v > abs(j - i):
               to_return.add(j)
               break

           dual_lambda_u = u_cum_dual[end] - u_cum_dual[start - 1]
           dual_lambda_v = v_cum_dual[end] - v_cum_dual[start - 1]

           if dual_lambda_u + dual_lambda_v - _sage_const_1 > abs(j - i):
               to_return.add(j)
               break

    return to_return

# =============================================================================
# STATISTICS FUNCTIONS (Unchanged)
# =============================================================================

def compute_expectation(frequencies):
    total_count = sum(frequencies.values())
    if total_count == 0:
        return 0.0
    weighted_sum = sum(value * count for value, count in frequencies.items())
    return weighted_sum / total_count

def compute_statistics(frequencies):
    if not frequencies: return {}
    total_count = sum(frequencies.values())
    expectation = compute_expectation(frequencies)
    mean_sq = sum(value**2 * count for value, count in frequencies.items()) / total_count
    variance = mean_sq - expectation**2
    std_dev = variance ** 0.5
    mode_value = max(frequencies.items(), key=lambda x: x[1])[0]
    mode_count = frequencies[mode_value]
    min_value = min(frequencies.keys())
    max_value = max(frequencies.keys())
    return {
        'expectation': expectation, 'variance': variance, 'std_dev': std_dev,
        'mode': mode_value, 'mode_count': mode_count, 'min': min_value,
        'max': max_value, 'total_count': total_count
    }

# =============================================================================
# MAIN SCRIPT (Modified for Grassmannian Permutations)
# Calculates the 'integer support' (Forward Stability Number, FS(u,v))
# for pairs of GRASSMANNIAN permutations in Gr(n) x Gr(n) and generates
# a bar chart of the frequencies using multiprocessing.
# =============================================================================

def generate_stability_chart(n, num_samples=None):
    # Generate all Grassmannian permutations in S_n (with caching)
    print(f"Getting Grassmannian permutations in S_{n}...")
    grass_perms = get_grassmannian_perms(n, use_cache=True)
    num_grass = len(grass_perms)
    expected = 2**n - n
    print(f"Total: {num_grass} Grassmannian permutations (expected: {expected})")

    # Verify count
    if num_grass != expected:
        print(f"WARNING: Count mismatch! Got {num_grass}, expected {expected}")

    task_list = [] # This will hold all the (u, v) pairs to process

    # Determine if we're sampling or computing all pairs
    if num_samples is None or num_samples >= num_grass ** 2:
        # Compute all pairs
        use_sampling = False
        num_to_process = num_grass ** 2
        print(f"Generating stability data for ALL pairs in Gr({n}) x Gr({n})...")
        print(f"Total pairs to process: {num_to_process}")
        print("Building task list...")
        task_list = [(u, v) for u in grass_perms for v in grass_perms]

    else:
        # Random sampling
        use_sampling = True
        num_to_process = num_samples
        print(f"Generating stability data using RANDOM SAMPLING from Gr({n}) x Gr({n})...")
        print(f"Sampling {num_to_process} random pairs")

        # Pre-generate pool of permutations for faster sampling
        print(f"Pre-generating pool of {min(POOL_SIZE, num_to_process * 2)} Grassmannian permutations...")
        pool_size = min(POOL_SIZE, num_to_process * 2, num_grass)

        # If pool size is larger than number of Grassmannian perms, just use all of them
        if pool_size >= num_grass:
            perm_pool = grass_perms
            pool_size = num_grass
        else:
            # Random sample from Grassmannian perms
            perm_pool = random.sample(grass_perms, pool_size)
        print(f"Pool generation complete with {pool_size} permutations.")

        # Generate random indices using NumPy (ultra-fast)
        print(f"Generating {num_to_process} random index pairs...")
        indices_u = np.random.randint(0, pool_size, num_to_process)
        indices_v = np.random.randint(0, pool_size, num_to_process)

        print("Building task list from random pairs...")
        # Use list comprehension for fast creation of the task list
        task_list = [(perm_pool[indices_u[i]], perm_pool[indices_v[i]]) for i in range(num_to_process)]
        print(f"Task list of {len(task_list)} pairs built.")

    # --- MULTIPROCESSING EXECUTION ---
    num_cores = cpu_count()
    print(f"\nStarting parallel processing with {num_cores} cores...")
    start_time = time.time()

    # Create the pool and distribute the work
    # Using imap with tqdm for progress tracking
    with Pool(processes=num_cores) as pool:
        # 'results' will be a list of integers: [4, 5, 4, 6, 7, 5, ...]
        chunksize = max(1, len(task_list) // (num_cores * 4))
        results = list(tqdm(
            pool.imap(process_pair, task_list, chunksize=chunksize),
            total=len(task_list),
            desc="Processing pairs",
            unit="pair"
        ))

    end_time = time.time()
    print(f"\nParallel processing finished in {end_time - start_time:.4f} seconds.")

    # --- AGGREGATE RESULTS ---
    print("Aggregating results...")
    # This is now incredibly simple: just count the results
    support_frequencies = Counter(results)

    print("\nData generation complete.")
    print(f"Frequencies: {support_frequencies}")

    if use_sampling:
        print(f"Note: These are estimated frequencies from {num_to_process} random samples")

    # --- STATISTICS AND CHARTING (Unchanged) ---
    print("\n" + "="*60)
    print("STATISTICAL ANALYSIS")
    print("="*60)

    stats = compute_statistics(support_frequencies)

    print(f"Total pairs analyzed: {stats['total_count']}")
    print(f"\nExpected Value E[FS(u,v)]: {stats['expectation']:.4f}")
    print(f"Standard Deviation: {stats['std_dev']:.4f}")
    print(f"Variance: {stats['variance']:.4f}")
    print(f"\nMode (most common): {stats['mode']} (appears {stats['mode_count']} times)")
    print(f"Minimum value: {stats['min']}")
    print(f"Maximum value: {stats['max']}")
    print("="*60)

    sorted_items = sorted(support_frequencies.items())
    x_values = [item[0] for item in sorted_items]
    y_values = [item[1] for item in sorted_items]

    print("Generating bar chart...")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x_values, y_values, color='steelblue', edgecolor='black')
    ax.set_xlabel('Integer Support (FS(u,v))', fontsize=12)
    ax.set_ylabel('Number of Pairs (Frequency)', fontsize=12)

    expectation = stats['expectation']
    std_dev = stats['std_dev']
    variance = stats['variance']
    title = f'Forward Stability Number for Grassmannian Permutations in S_{n}\n'
    title += f'E[FS(u,v)] = {expectation:.4f}, σ = {std_dev:.4f}, σ² = {variance:.4f}'
    ax.set_title(title, fontsize=14)

    ax.set_xticks(range(min(x_values), max(x_values) + 1)) # Ensure all integer ticks
    ax.grid(axis='y', alpha=0.3)

    ax.axvline(x=expectation, color='red', linestyle='--', linewidth=2,
               label=f'Expected Value = {expectation:.4f}', alpha=0.7)
    ax.legend()

    if use_sampling:
        chart_filename = f"grassmannian_s{n}_stability_sampled_{num_to_process}.png"
    else:
        chart_filename = f"grassmannian_s{n}_stability.png"

    plt.tight_layout()
    plt.savefig(chart_filename, dpi=150)
    plt.close()

    print(f"\nChart saved as '{chart_filename}'")
    print("Script finished.")

    return support_frequencies, stats

# --- Execute the function ---
if __name__ == "__main__":
    # This check is essential for multiprocessing to work!
    generate_stability_chart(n=N, num_samples=NUM_RANDOM_SAMPLES)
