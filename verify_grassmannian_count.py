#!/usr/bin/env python3
"""
Verify the combinatorial logic for Grassmannian permutation generation
(without needing SageMath)
"""
from itertools import combinations

def count_grassmannian_direct(n):
    """
    Count Grassmannian permutations using direct generation logic.
    Returns the count without creating actual permutation objects.
    """
    count = 1  # Identity permutation (no descent)

    # For each descent position k
    for k in range(1, n):
        # Choose which k elements appear in first k positions
        for first_subset in combinations(range(1, n+1), k):
            first_part = list(first_subset)
            second_part = sorted(set(range(1, n+1)) - set(first_part))

            # Check if this creates a descent: max(first) > min(second)
            if first_part[-1] > second_part[0]:
                count += 1

    return count

print("="*60)
print("GRASSMANNIAN PERMUTATION COUNT VERIFICATION")
print("="*60)
print("\nVerifying formula: |Gr(n)| = 2^n - n")
print()

for n in range(1, 15):
    actual = count_grassmannian_direct(n)
    expected = 2**n - n
    status = "✓" if actual == expected else "✗"
    print(f"n={n:2d}: actual={actual:5d}, expected={expected:5d} {status}")

print("\n" + "="*60)
print("All counts verified!" if all(
    count_grassmannian_direct(n) == 2**n - n for n in range(1, 15)
) else "ERROR: Some counts don't match!")
print("="*60)
