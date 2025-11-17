from fractions import Fraction
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
    Compute the theoretical formula: 2^(j-1) - sum_{k=0}^{j-1} C(n,k)
    
    Args:
        n: Size of permutation
        j: Position (1-indexed)
    
    Returns:
        The value of the formula
    """
    # Compute sum of binomial coefficients
    binom_sum = sum(binomial(n, k) for k in range(j))
    
    # Compute 2^(j-1) - binom_sum
    result = 2**(j-1) - binom_sum
    
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
    formula_value = int(compute_formula(n, j))
    total_expected = 2**n - n
    normalized_count = Fraction(count, total_expected)
    normalized_formula = Fraction(formula_value, total_expected)
    normalized_difference = normalized_count - normalized_formula
    matches = (normalized_difference == 0)
    
    if verbose:
        print(f"=" * 70)
        print(f"Analysis for n = {n}, position j = {j}")
        print(f"=" * 70)
        print(f"Total Grassmannian permutations (enumerated): {total_grass}")
        print(f"Expected (2^n - n): {total_expected}")
        if total_grass != total_expected:
            print(f"WARNING: Enumerated total differs from 2^n - n by {total_grass - total_expected}.")
        print(f"\nPermutations with LTR max at position {j}: {count}")
        print(f"Probability P(d_{j} = 0) = {count}/{total_expected} = {float(normalized_count):.4f}")
        print(f"Compare to 1/{j} = {1/j:.4f} (full symmetric group)")
        
        # Formula verification
        print(f"\n{'=' * 70}")
        print(f"FORMULA VERIFICATION")
        print(f"{'=' * 70}")
        print(f"Formula: 2^(j-1) - sum_{{k=0}}^{{j-1}} C(n,k)")
        print(f"       = 2^({j}-1) - sum_{{k=0}}^{{{j-1}}} C({n},k)")
        
        # Show the binomial sum computation
        binom_sum = sum(binomial(n, k) for k in range(j))
        print(f"       = {2**(j-1)} - {binom_sum}")
        print(f"       = {formula_value}")
        print(f"\nActual count: {count}")
        print(f"Formula value: {formula_value}")
        print(f"\nNormalized count  = {normalized_count} (~{float(normalized_count):.6f})")
        print(f"Normalized formula = {normalized_formula} (~{float(normalized_formula):.6f})")
        
        if matches:
            print(f"✓ MATCH after dividing by (2^n - n).")
        else:
            print(f"✗ NO MATCH after dividing by (2^n - n).")
            print(f"  Normalized difference = {normalized_difference} "
                  f"(~{float(normalized_difference):+.6f})")
        
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
        'probability': float(normalized_count),
        'formula_value': formula_value,
        'normalized_count': normalized_count,
        'normalized_formula': normalized_formula,
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
        print(f"\n" + "=" * 135)
        print(f"Summary for n = {n} (Total Grassmannian permutations: {results[1]['total_grassmannian']})")
        print(f"=" * 135)
        print(f"{'Pos j':^6} | {'Count':^6} | {'Formula':^8} | {'Count/(2^n-n)':^18} | "
              f"{'Formula/(2^n-n)':^20} | {'Match?':^7} | {'1/j':^10} | {'Diff (count-1/j)':^18}")
        print(f"{'-' * 135}")
        
        all_match = True
        for j in range(1, n+1):
            count = results[j]['count_with_ltr_max']
            formula_val = results[j]['formula_value']
            matches = results[j]['matches_formula']
            prob = results[j]['probability']
            count_norm = results[j]['normalized_count']
            formula_norm = results[j]['normalized_formula']
            expected = 1/j
            diff = prob - expected
            
            match_str = "✓" if matches else "✗"
            if not matches:
                all_match = False
            
            print(f"{j:^6d} | {count:^6d} | {formula_val:^8d} | {float(count_norm):^18.6f} | "
                  f"{float(formula_norm):^20.6f} | {match_str:^7s} | "
                  f"{expected:^10.4f} | {diff:^+18.6f}")
        
        print(f"=" * 135)
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
        print(f"Formula: 2^(j-1) - sum_{{k=0}}^{{j-1}} C(n,k)")
        print(f"Comparison performed on normalized values value / (2^n - n)")
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
                        count_norm = results[j]['normalized_count']
                        formula_norm = results[j]['normalized_formula']
                        print(f"    Position {j}:")
                        print(f"      Count            = {count}")
                        print(f"      Formula value    = {formula_val}")
                        print(f"      Count/(2^n - n)  = {count_norm} (~{float(count_norm):.6f})")
                        print(f"      Formula/(2^n - n)= {formula_norm} (~{float(formula_norm):.6f})")
    
    return all_results

def detailed_formula_breakdown(n, j):
    """
    Show detailed breakdown of the formula computation.
    
    Args:
        n: Size of permutation
        j: Position
    """
    print(f"\nDetailed Formula Breakdown for n={n}, j={j}")
    print(f"=" * 60)
    print(f"Formula: 2^(j-1) - sum_{{k=0}}^{{j-1}} C(n,k)")
    print(f"\nStep 1: Compute 2^(j-1) = 2^{j-1} = {2**(j-1)}")
    print(f"\nStep 2: Compute sum_{{k=0}}^{{{j-1}}} C({n},k)")
    
    total_sum = 0
    for k in range(j):
        binom_val = binomial(n, k)
        total_sum += binom_val
        print(f"        C({n},{k}) = {binom_val}")
    
    print(f"        Sum = {total_sum}")
    formula_value = 2**(j-1) - total_sum
    print(f"\nStep 3: Final result = {2**(j-1)} - {total_sum} = {formula_value}")
    
    # Verify with actual count
    perms = find_grassmannian_with_ltr_max(n, j)
    count = len(perms)
    print(f"\nActual count from enumeration: {count}")
    denom = 2**n - n
    normalized_count = Fraction(count, denom)
    normalized_formula = Fraction(formula_value, denom)
    print(f"Raw match (counts)       : {count == formula_value}")
    print(f"Normalized count         : {normalized_count} (~{float(normalized_count):.6f})")
    print(f"Normalized formula value : {normalized_formula} (~{float(normalized_formula):.6f})")
    print(f"Match after normalization: {normalized_count == normalized_formula}")

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
print("EXAMPLE 3: Formula verification for n=3,4,5,6,7")
verify_formula_range(range(3, 8))

print("\n" * 2)

# Example 4: Detailed formula breakdown
print("EXAMPLE 4: Detailed formula breakdown")
detailed_formula_breakdown(5, 3)
detailed_formula_breakdown(6, 4)

print("\n" * 2)

# Example 5: Quick formula check function
def quick_check(n, j):
    """Quick check if count matches formula for given n and j"""
    perms = find_grassmannian_with_ltr_max(n, j)
    count = len(perms)
    formula_val = int(compute_formula(n, j))
    denom = 2**n - n
    normalized_count = Fraction(count, denom)
    normalized_formula = Fraction(formula_val, denom)
    matches = normalized_count == normalized_formula
    
    print(f"n={n}, j={j}: Count={count}, Formula={formula_val}")
    print(f"           Count/(2^n-n)={normalized_count} (~{float(normalized_count):.6f}), "
          f"Formula/(2^n-n)={normalized_formula} (~{float(normalized_formula):.6f}), "
          f"Match={matches}")
    return matches

# Test a few cases
print("EXAMPLE 5: Quick checks")
for n in [4, 5, 6]:
    for j in [1, 2, 3]:
        if j <= n:
            quick_check(n, j)
