# Import necessary libraries
from sage.all import *
from collections import Counter
import random
import numpy as np

# --- Define constants ---
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)

# --- Configuration ---
# Change this value to compute for different S_n
N = 60

# Random sampling configuration
# Set to None to compute ALL pairs, or set to a number to randomly sample that many pairs
NUM_RANDOM_SAMPLES = 100000  # e.g., 100000 random pairs instead of all pairs

# Pool size for pre-generated permutations (for faster sampling)
POOL_SIZE = 7000  # Medium pool size for balance of speed and diversity

# Random seed for reproducibility (set to None for different results each run)
RANDOM_SEED = None

# --- Helper functions from your provided code ---

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
    
    Returns two lists:
    1. cum_theta: cum_theta[i] = number of nonzero Lehmer entries from 1 to i
    2. cum_dual_theta: cum_dual_theta[i] = number of nonzero Dual Lehmer entries from 1 to i
    """
    n = len(w)
    
    # 1. Get Sage's optimized Lehmer and Dual Lehmer codes
    # These are 0-indexed lists of length n
    lehmer_code = w.to_lehmer_code()
    dual_lehmer_code = w.to_lehmer_cocode()

    # 2. Create cumulative sum arrays, 1-indexed for easy lookup
    # We make them of size max_k + 1
    cum_theta = [0] * (max_k + 1)
    cum_dual_theta = [0] * (max_k + 1)
    
    for k in range(1, n + 1):
        # Fill cumulative arrays
        cum_theta[k] = cum_theta[k - 1]
        cum_dual_theta[k] = cum_dual_theta[k - 1]
        
        if lehmer_code[k - 1] > 0:
            cum_theta[k] += 1
            
        if dual_lehmer_code[k - 1] > 0:
            cum_dual_theta[k] += 1
            
    # Propagate the last value for indices k > n
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
    n = len(u) # Assumes len(u) == len(v)
    
    # Set the maximum index to check
    max_n = 2 * n + _sage_const_1
    
    # 1. Precompute all lambda values for u and v
    u_cum_theta, u_cum_dual = precompute_lambda_lookups(u, max_n)
    v_cum_theta, v_cum_dual = precompute_lambda_lookups(v, max_n)

    # 2. Add base supports
    to_return.update(integer_support(u))
    to_return.update(integer_support(v))

    # 3. Run the O(n^2) loop
    for j in range(_sage_const_1, max_n):
        if j in to_return:
            continue
            
        for i in range(_sage_const_1, max_n):
           # Calculate ordinary_lambda(w, j, i) in O(1)
           # This is sum(theta[k] for k from i to j).
           # Note: Your original code had i and j swapped here.
           # ordinary_lambda(u, j, i) implies j is start, i is end.
           start, end = min(i, j), max(i, j)
           ord_lambda_u = u_cum_theta[end] - u_cum_theta[start - 1]
           ord_lambda_v = v_cum_theta[end] - v_cum_theta[start - 1]

           if ord_lambda_u + ord_lambda_v > abs(j - i):
               to_return.add(j)
               break
               
           # Calculate dual_lambda(w, i, j) in O(1)
           # This is sum(dual_theta[k] for k from i to j).
           start, end = min(i, j), max(i, j)
           dual_lambda_u = u_cum_dual[end] - u_cum_dual[start - 1]
           dual_lambda_v = v_cum_dual[end] - v_cum_dual[start - 1]
           
           if dual_lambda_u + dual_lambda_v - _sage_const_1 > abs(j - i):
               to_return.add(j)
               break
               
    return to_return

def compute_expectation(frequencies):
    """
    Computes the expected value (average) of the forward stability number
    from a frequency distribution.
    
    Args:
        frequencies: A Counter or dict with {stability_value: count}
    
    Returns:
        float: The expected value E[FS(u,v)]
    """
    total_count = sum(frequencies.values())
    if total_count == 0:
        return 0.0
    
    weighted_sum = sum(value * count for value, count in frequencies.items())
    expectation = weighted_sum / total_count
    
    return expectation

def compute_statistics(frequencies):
    """
    Computes various statistics from the frequency distribution.
    
    Args:
        frequencies: A Counter or dict with {stability_value: count}
    
    Returns:
        dict: Dictionary containing statistical measures
    """
    if not frequencies:
        return {}
    
    total_count = sum(frequencies.values())
    expectation = compute_expectation(frequencies)
    
    # Compute variance and standard deviation
    mean_sq = sum(value**2 * count for value, count in frequencies.items()) / total_count
    variance = mean_sq - expectation**2
    std_dev = variance ** 0.5
    
    # Find mode (most common value)
    mode_value = max(frequencies.items(), key=lambda x: x[1])[0]
    mode_count = frequencies[mode_value]
    
    # Find min and max
    min_value = min(frequencies.keys())
    max_value = max(frequencies.keys())
    
    stats = {
        'expectation': expectation,
        'variance': variance,
        'std_dev': std_dev,
        'mode': mode_value,
        'mode_count': mode_count,
        'min': min_value,
        'max': max_value,
        'total_count': total_count
    }
    
    return stats

# --- Main script ---

def generate_stability_chart(n, num_samples=None, random_seed=None):
    """
    Calculates the 'integer support' (Forward Stability Number, FS(u,v))
    for pairs of permutations in S_n x S_n and generates a bar chart
    of the frequencies.
    
    Args:
        n: The size of the symmetric group S_n
        num_samples: If None, compute all pairs. If an integer, randomly sample that many pairs.
        random_seed: Seed for random number generator (for reproducibility)
    """
    # Set random seed if provided
    if random_seed is not None:
        random.seed(random_seed)
    
    # Use Counter to store frequencies
    support_frequencies = Counter()
    
    # Create Permutations object once for efficiency
    perms_obj = Permutations(n)
    
    # Determine if we're sampling or computing all pairs
    if num_samples is None or num_samples >= factorial(n) ** 2:
        # Compute all pairs
        use_sampling = False
        perms = list(perms_obj)
        num_to_process = len(perms) ** 2
        print(f"Generating stability data for ALL pairs in S_{n} x S_{n}...")
        print(f"Total pairs to process: {num_to_process}")
    else:
        # Random sampling
        use_sampling = True
        num_to_process = num_samples
        print(f"Generating stability data using RANDOM SAMPLING from S_{n} x S_{n}...")
        print(f"Total possible pairs: {factorial(n) ** 2}")
        print(f"Sampling {num_to_process} random pairs")

    if use_sampling:
        # Pre-generate pool of permutations for faster sampling
        print(f"Pre-generating pool of {min(POOL_SIZE, num_to_process * 2)} permutations...")
        pool_size = min(POOL_SIZE, num_to_process * 2)  # Don't generate more than needed
        perm_pool = [perms_obj.random_element() for _ in range(pool_size)]
        print(f"Pool generation complete.")
        
        # Generate random pairs using NumPy (ultra-fast batch generation)
        print(f"Generating {num_to_process} random index pairs...")
        final_indices = np.empty((num_to_process, 2), dtype=np.intp)
        count = 0
        
        while count < num_to_process:
            # How many more do we need?
            needed = num_to_process - count
            
            # Generate a new batch of random indices
            i = np.random.randint(0, pool_size, needed)
            j = np.random.randint(0, pool_size, needed)
            
            # Create a filter (mask) for i != j (avoid self-pairs if desired)
            # Note: For this application, i == j is actually valid (same permutation paired with itself)
            # but we'll keep the filter as optional
            # If you want to allow i == j, simply use: mask = np.ones(needed, dtype=bool)
            mask = (i != j)
            
            # Get the valid indices
            valid_i = i[mask]
            valid_j = j[mask]
            
            # Stack them into pairs
            valid_pairs = np.stack((valid_i, valid_j), axis=1)
            
            # How many valid pairs did we get?
            num_valid = len(valid_pairs)
            if num_valid == 0:  # In case of bad luck
                continue
            
            # How many can we add to our final array?
            num_to_add = min(needed, num_valid)
            
            # Add them
            final_indices[count : count + num_to_add] = valid_pairs[:num_to_add]
            count += num_to_add
        
        print(f"Index generation complete. Processing {num_to_process} pairs...")
        
        # Now process all pairs using the pre-generated indices
        for pairs_processed in range(num_to_process):
            # Get the pair of permutations using list indexing (keeps Sage objects intact)
            idx_u, idx_v = final_indices[pairs_processed]
            u = perm_pool[idx_u]
            v = perm_pool[idx_v]
            
            # Get the support set from the formula
            support_set = integer_support_from_formula(u, v)
            
            # Find the largest number in the set
            integer_support_val = max(support_set) if support_set else 0
            
            # Increment the count for this value
            support_frequencies[integer_support_val] += 1
            
            # Print progress (reduced frequency for speed)
            if (pairs_processed + 1) % max(1, num_to_process // 10) == 0:
                print(f"Processed {pairs_processed + 1} / {num_to_process} pairs...")
    else:
        # Compute all pairs (original approach)
        for i, u in enumerate(perms):
            for v in perms:
                # Get the support set from the formula
                support_set = integer_support_from_formula(u, v)
                
                # Find the largest number in the set (the Forward Stability Number, FS(u,v))
                # If the set is empty, the stability number is 0
                integer_support_val = max(support_set) if support_set else 0
                
                # Increment the count for this value
                support_frequencies[integer_support_val] += 1
                
            # Optional: print progress
            if (i + 1) % 6 == 0:
                print(f"Processed { (i + 1) * len(perms) } / {num_to_process} pairs...")

    print("\nData generation complete.")
    print(f"Frequencies: {support_frequencies}")
    
    if use_sampling:
        print(f"Note: These are estimated frequencies from {num_to_process} random samples")
        print(f"      Sampled uniformly from pool of {pool_size} pre-generated permutations")
    
    # Compute and display statistics
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
    
    # Prepare data for the bar chart
    # Sort by the integer support value (the key)
    sorted_items = sorted(support_frequencies.items())
    
    # Extract x-values (integer support) and y-values (frequencies)
    x_values = [item[0] for item in sorted_items]
    y_values = [item[1] for item in sorted_items]
    
    print("Generating bar chart...")
    
    # Use matplotlib directly through Sage for better control
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x_values, y_values, color='steelblue', edgecolor='black')
    ax.set_xlabel('Integer Support (FS(u,v))', fontsize=12)
    ax.set_ylabel('Number of Pairs (Frequency)', fontsize=12)
    
    # Use already computed expectation from stats
    expectation = stats['expectation']
    title = f'Forward Stability Number for S_{n}\n'
    title += f'E[FS(u,v)] = {expectation:.4f}'
    ax.set_title(title, fontsize=14)
    
    ax.set_xticks(x_values)
    ax.grid(axis='y', alpha=0.3)
    
    # Add a vertical line at the expected value
    ax.axvline(x=expectation, color='red', linestyle='--', linewidth=2, 
               label=f'Expected Value = {expectation:.4f}', alpha=0.7)
    ax.legend()
    
    # Save the chart to a file
    if use_sampling:
        chart_filename = f"s{n}_stability_frequencies_sampled_{num_to_process}.png"
    else:
        chart_filename = f"s{n}_stability_frequencies.png"
    
    plt.tight_layout()
    plt.savefig(chart_filename, dpi=150)
    plt.close()
    
    print(f"\nChart saved as '{chart_filename}'")
    print("Script finished.")
    
    return support_frequencies, stats

# --- Execute the function ---
if __name__ == "__main__":
    generate_stability_chart(n=N, num_samples=NUM_RANDOM_SAMPLES, random_seed=RANDOM_SEED)
