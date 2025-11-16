from sage.all import *

def is_left_to_right_max(perm, j):
    """
    Check if position j is a left-to-right maximum in permutation perm.
    
    Args:
        perm: A permutation (1-indexed in Sage)
        j: Position to check (1-indexed)
    
    Returns:
        True if w(j) is larger than all w(k) for k < j
    """
    if j == 1:
        return True
    
    # Check if all elements before position j are smaller than perm[j-1]
    value_at_j = perm[j-1]  # Sage permutations are 0-indexed internally
    
    for k in range(j-1):
        if perm[k] >= value_at_j:
            return False
    
    return True

def compute_formula(n, j):
    """
    Compute the theoretical formula: sum_{k=0}^{n-j} C(n,k) + 2^(j-1) - n
    
    Args:
        n: Size of permutation
        j: Position (1-indexed)
    
    Returns:
        The value of the formula
    """
    # Compute sum of binomial coefficients
    binom_sum = sum(binomial(n, k) for k in range(n - j + 1))
    
    # Add 2^(j-1) and subtract n
    result = binom_sum + 2^(j-1) - n
    
    return result

def find_grassmannian_with_ltr_max(n, j):
    """
    Find all Grassmannian permutations where position j is a left-to-right maximum.
    
    Args:
        n: Size of permutation
        j: Position to check (1-indexed, must satisfy 1 <= j <= n)
    
    Returns:
        List of Grassmannian permutations with LTR max at position j
    """
    if not (1 <= j <= n):
        raise ValueError(f"j must be between 1 and {n}")
    
    # Generate all permutations and filter for Grassmannian
    P = Permutations(n)
    grassmannian_perms = [p for p in P if p.number_of_descents() <= 1]
    
    # Filter for those with LTR max at position j
    result = [p for p in grassmannian_perms if is_left_to_right_max(p, j)]
    
    return result

def analyze_grassmannian_ltr_max(n, j, verbose=True):
    """
    Analyze Grassmannian permutations with LTR max at position j.
    Includes verification against the theoretical formula.
    
    Args:
        n: Size of permutation
        j: Position to check
        verbose: If True, print detailed information
    
    Returns:
        Dictionary with analysis results
    """
    # Get all Grassmannian permutations
    P = Permutations(n)
    all_grass = [p for p in P if p.number_of_descents() <= 1]
    total_grass = len(all_grass)
    
    # Get those with LTR max at position j
    ltr_max_at_j = find_grassmannian_with_ltr_max(n, j)
    count = len(ltr_max_at_j)
    
    # Compute theoretical formula
    formula_value = compute_formula(n, j)
    matches = (count == formula_value)
    
    if verbose:
        print(f"=" * 70)
        print(f"Analysis for n = {n}, position j = {j}")
        print(f"=" * 70)
        print(f"Total Grassmannian permutations: {total_grass}")
        print(f"Expected (2^n - n): {2**n - n}")
        print(f"\nPermutations with LTR max at position {j}: {count}")
        print(f"Probability P(d_{j} = 0) = {count}/{total_grass} = {count/total_grass:.4f}")
        print(f"Compare to 1/{j} = {1/j:.4f} (full symmetric group)")
        
        # Formula verification
        print(f"\n{'=' * 70}")
        print(f"FORMULA VERIFICATION")
        print(f"{'=' * 70}")
        print(f"Formula: sum_{{k=0}}^{{n-j}} C(n,k) + 2^(j-1) - n")
        print(f"       = sum_{{k=0}}^{{{n-j}}} C({n},k) + 2^{j-1} - {n}")
        
        # Show the binomial sum computation
        binom_sum = sum(binomial(n, k) for k in range(n - j + 1))
        print(f"       = {binom_sum} + {2^(j-1)} - {n}")
        print(f"       = {formula_value}")
        print(f"\nActual count: {count}")
        print(f"Formula value: {formula_value}")
        
        if matches:
            print(f"✓ MATCH! The formula is correct.")
        else:
            print(f"✗ NO MATCH! Difference: {count - formula_value}")
        
        print(f"\n{'-' * 70}")
        print(f"All Grassmannian permutations with LTR max at position {j}:")
        print(f"{'-' * 70}")
        
        for i, p in enumerate(ltr_max_at_j, 1):
            descents = p.descents()
            inversions = p.number_of_inversions()
            
            # Show the permutation in one-line notation
            perm_str = str(list(p))
            
            # Mark position j
            desc_str = f"at position {descents[0]}" if descents else "none"
            
            print(f"{i:3d}. {perm_str:20s} | Descent: {desc_str:20s} | "
                  f"Inversions: {inversions:2d} | w({j}) = {p[j-1]}")
    
    return {
        'total_grassmannian': total_grass,
        'count_with_ltr_max': count,
        'probability': count / total_grass,
        'formula_value': formula_value,
        'matches_formula': matches,
        'permutations': ltr_max_at_j
    }

def analyze_all_positions(n, verbose=True):
    """
    Analyze LTR maxima at all positions for Grassmannian permutations.
    Includes formula verification for each position.
    
    Args:
        n: Size of permutation
        verbose: If True, print summary table
    
    Returns:
        Dictionary mapping position j to analysis results
    """
    results = {}
    
    for j in range(1, n+1):
        results[j] = analyze_grassmannian_ltr_max(n, j, verbose=False)
    
    if verbose:
        print(f"\n" + "=" * 95)
        print(f"Summary for n = {n} (Total Grassmannian permutations: {results[1]['total_grassmannian']})")
        print(f"=" * 95)
        print(f"{'Pos j':^6} | {'Count':^6} | {'Formula':^8} | {'Match?':^7} | "
              f"{'P(d_j=0)':^10} | {'1/j':^10} | {'Difference':^10}")
        print(f"{'-' * 95}")
        
        all_match = True
        for j in range(1, n+1):
            count = results[j]['count_with_ltr_max']
            formula_val = results[j]['formula_value']
            matches = results[j]['matches_formula']
            total = results[j]['total_grassmannian']
            prob = results[j]['probability']
            expected = 1/j
            diff = prob - expected
            
            match_str = "✓" if matches else "✗"
            if not matches:
                all_match = False
            
            print(f"{j:^6d} | {count:^6d} | {formula_val:^8d} | {match_str:^7s} | "
                  f"{prob:^10.4f} | {expected:^10.4f} | {diff:^+10.4f}")
        
        print(f"=" * 95)
        if all_match:
            print(f"✓ All positions match the formula!")
        else:
            print(f"✗ Some positions do not match the formula.")
    
    return results

def verify_formula_range(n_values, verbose=True):
    """
    Verify the formula across multiple values of n.
    
    Args:
        n_values: List or range of n values to test
        verbose: If True, print results
    
    Returns:
        Dictionary with verification results
    """
    all_results = {}
    
    for n in n_values:
        all_results[n] = analyze_all_positions(n, verbose=False)
    
    if verbose:
        print(f"\n" + "=" * 80)
        print(f"FORMULA VERIFICATION ACROSS MULTIPLE VALUES OF n")
        print(f"=" * 80)
        
        for n in n_values:
            print(f"\nn = {n}:")
            results = all_results[n]
            
            # Check if all positions match
            all_match = all(results[j]['matches_formula'] for j in range(1, n+1))
            
            if all_match:
                print(f"  ✓ Formula verified for ALL positions (j = 1, 2, ..., {n})")
            else:
                print(f"  ✗ Formula does NOT match for some positions:")
                for j in range(1, n+1):
                    if not results[j]['matches_formula']:
                        count = results[j]['count_with_ltr_max']
                        formula_val = results[j]['formula_value']
                        print(f"    Position {j}: Count = {count}, Formula = {formula_val}")
    
    return all_results

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Example 1: Detailed analysis for a single case
print("EXAMPLE 1: Detailed analysis for n=5, j=3")
analyze_grassmannian_ltr_max(5, 3)

print("\n" * 2)

# Example 2: Summary table for all positions at n=5
print("EXAMPLE 2: Summary for n=5, all positions")
analyze_all_positions(5)

print("\n" * 2)

# Example 3: Verify formula across multiple n values
print("EXAMPLE 3: Formula verification for n=3,4,5,6")
verify_formula_range(range(3, 7))

print("\n" * 2)

# Example 4: Quick formula check function
def quick_check(n, j):
    """Quick check if count matches formula for given n and j"""
    perms = find_grassmannian_with_ltr_max(n, j)
    count = len(perms)
    formula_val = compute_formula(n, j)
    
    print(f"n={n}, j={j}: Count={count}, Formula={formula_val}, Match={count==formula_val}")
    return count == formula_val

# Test a few cases
print("EXAMPLE 4: Quick checks")
for n in [4, 5, 6]:
    for j in [1, 2, 3]:
        if j <= n:
            quick_check(n, j)
```

This enhanced version includes:

1. **`compute_formula(n, j)`**: Computes the theoretical formula value
2. **Formula verification** in the detailed analysis showing:
   - The formula with actual values substituted
   - Step-by-step computation
   - Whether it matches the actual count
3. **Enhanced summary table** showing both count and formula value with a match indicator (✓ or ✗)
4. **`verify_formula_range()`**: Tests the formula across multiple n values
5. **`quick_check()`**: A simple function for rapid verification

**Sample Output:**
```
======================================================================
FORMULA VERIFICATION
======================================================================
Formula: sum_{k=0}^{n-j} C(n,k) + 2^(j-1) - n
       = sum_{k=0}^{2} C(5,k) + 2^2 - 5
       = 16 + 4 - 5
       = 15

Actual count: 15
Formula value: 15
✓ MATCH! The formula is correct.
