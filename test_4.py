# Import necessary libraries
from sage.all import *
from collections import Counter
import random

# --- Define constants ---
_sage_const_1 = Integer(1); _sage_const_0 = Integer(0);

# --- Configuration ---
# Change this value to compute for different S_n
N = 50

# Random sampling configuration
# Set to None to compute ALL pairs, or set to a number to randomly sample that many pairs
NUM_RANDOM_SAMPLES = 5000  # e.g., 10000 random pairs instead of all pairs

# Random seed for reproducibility (set to None for different results each run)
RANDOM_SEED = None

# --- Helper functions from your provided code ---

def ordinary_lambda(w, i, j):
    """
    Calculates the number of nonzero Lehmer code entries
    for permutation w between indices i and j (inclusive).
    """
    lehmer_theta = []
    for k in range(i, j + _sage_const_1):
        # Ensure we don't index out of bounds for permutations
        if k > len(w):
            continue
        to_append = _sage_const_0
        for l in range(k + _sage_const_1, len(w) + _sage_const_1):
            if w[l - _sage_const_1] < w[k - _sage_const_1]:
                to_append = _sage_const_1
                break # Found one, no need to check further for this k
        lehmer_theta.append(to_append)
    return sum(lehmer_theta)

def dual_lambda(w, i, j):
    """
    Calculates the number of nonzero dual Lehmer code entries
    for permutation w between indices i and j (inclusive).
    """
    dual_lehmer_theta = []
    for k in range(i, min(len(w) + _sage_const_1, j + _sage_const_1)):
        to_append = _sage_const_0
        for l in range(_sage_const_1, k):
            if w[l - _sage_const_1] > w[k - _sage_const_1]:
                to_append = _sage_const_1
                break # Found one, no need to check further for this k
        dual_lehmer_theta.append(to_append)
    return sum(dual_lehmer_theta)

def integer_support(w):
    """
    Finds the set of indices i for which w is not in S_i.
    This corresponds to {i | FS(w) > i}.
    """
    to_return = set()
    for i in range(_sage_const_1, len(w)):
        for j in range(_sage_const_1, i + _sage_const_1):
            if w[j - _sage_const_1] > i:
                to_return.add(i)
                break # Found one, no need to check further for this i
    return to_return

def integer_support_from_formula(u, v):
    """
    Calculates the set of unstable positions for the product S_u * S_v
    using the formulas from the paper.
    """
    to_return = set()
    # The maximum possible stability number for S_n x S_n is 2n-1.
    # Checking up to 2*len(u) is sufficient.
    max_n = 2 * len(u) + _sage_const_1
    
    # Add base supports
    to_return.update(integer_support(u))
    to_return.update(integer_support(v))

    for j in range(_sage_const_1, max_n):
        # Optimization: if j is already in the set, skip inner loop
        if j in to_return:
            continue
            
        for i in range(_sage_const_1, max_n):
           # Check ordinary lambda condition from your code
           if ordinary_lambda(u, j, i) + ordinary_lambda(v, j, i) > abs(j - i):
               to_return.add(j)
               break # j is added, move to next j
               
           # Check dual lambda condition from your code
           if dual_lambda(u, i, j) + dual_lambda(v, i, j) - _sage_const_1 > abs(j - i):
               to_return.add(j)
               break # j is added, move to next j
               
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
        # Random sampling approach - use Sage's random_element() method
        for pairs_processed in range(num_to_process):
            # Randomly select two permutations using Sage's method
            u = perms_obj.random_element()
            v = perms_obj.random_element()
            
            # Get the support set from the formula
            support_set = integer_support_from_formula(u, v)
            
            # Find the largest number in the set (the Forward Stability Number, FS(u,v))
            # If the set is empty, the stability number is 0
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
