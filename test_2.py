# Import necessary libraries
from sage.all import *
from collections import Counter

# --- Define constants ---
_sage_const_1 = Integer(1); _sage_const_0 = Integer(0);

# --- Configuration ---
# Change this value to compute for different S_n
N = 6

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

# --- Main script ---

def generate_stability_chart(n):
    """
    Calculates the 'integer support' (Forward Stability Number, FS(u,v))
    for all pairs of permutations in S_n x S_n and generates a bar chart
    of the frequencies.
    """
    print(f"Generating stability data for all pairs in S_{n} x S_{n}...")
    
    # Use Counter to store frequencies
    support_frequencies = Counter()
    
    # Get all permutations for S_n
    perms = list(Permutations(n))
    total_pairs = len(perms) ** 2
    
    print(f"Total pairs to process: {total_pairs} (for S_{n})")

    # Iterate over all possible pairs (u, v)
    for i, u in enumerate(perms):
        for v in perms:
            # Get the support set from the formula
            support_set = integer_support_from_formula(u, v)
            
            # Find the largest number in the set (the Forward Stability Number, FS(u,v))
            # If the set is empty, the stability number is 0
            if support_set:
                integer_support_val = max(support_set)
            else:
                integer_support_val = _sage_const_0
            
            # Increment the count for this value
            support_frequencies[integer_support_val] += 1
            
        # Optional: print progress
        if (i + 1) % 6 == 0:
            print(f"Processed { (i + 1) * len(perms) } / {total_pairs} pairs...")

    print("Data generation complete.")
    print(f"Frequencies: {support_frequencies}")
    
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
    ax.set_title(f'Forward Stability Number Frequencies for S_{n} x S_{n}', fontsize=14)
    ax.set_xticks(x_values)
    ax.grid(axis='y', alpha=0.3)
    
    # Save the chart to a file
    chart_filename = f"s{n}_stability_frequencies.png"
    plt.tight_layout()
    plt.savefig(chart_filename, dpi=150)
    plt.close()
    
    print(f"Chart saved as '{chart_filename}'")
    print("Script finished.")

# --- Execute the function ---
if __name__ == "__main__":
    generate_stability_chart(n=N)