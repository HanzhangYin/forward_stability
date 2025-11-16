#!/usr/bin/env python3
"""
Quick test to verify the performance of the optimized Grassmannian generation
"""
from sage.all import *
import time
from itertools import combinations

def grassmannian_perms_direct(n):
    """FAST direct generation using combinatorial structure."""
    P = Permutations(n)
    result = []

    # Identity permutation (no descent)
    result.append(P(list(range(1, n+1))))

    # For each descent position k
    for k in range(1, n):
        # Choose which k elements appear in first k positions
        for first_subset in combinations(range(1, n+1), k):
            first_part = list(first_subset)
            second_part = sorted(set(range(1, n+1)) - set(first_part))

            # Check if this creates a descent: max(first) > min(second)
            if first_part[-1] > second_part[0]:
                perm_list = first_part + second_part
                result.append(P(perm_list))

    return result

def grassmannian_perms_slow(n):
    """OLD slow method: filter all permutations."""
    P = Permutations(n)
    return [p for p in P if len(p.descents()) <= 1]

# Test for small n to compare
print("="*60)
print("PERFORMANCE COMPARISON: Direct Generation vs Filtering")
print("="*60)

for n in [5, 7, 9]:
    expected = 2**n - n
    print(f"\nTesting n={n} (expected {expected} permutations):")

    # Test direct generation (FAST)
    print("  Direct generation:", end=" ", flush=True)
    start = time.time()
    direct = grassmannian_perms_direct(n)
    direct_time = time.time() - start
    print(f"{len(direct)} perms in {direct_time:.4f}s")

    # Test filtering (SLOW)
    print("  Filtering method: ", end=" ", flush=True)
    start = time.time()
    filtered = grassmannian_perms_slow(n)
    filter_time = time.time() - start
    print(f"{len(filtered)} perms in {filter_time:.4f}s")

    # Speedup
    speedup = filter_time / direct_time
    print(f"  Speedup: {speedup:.1f}x faster!")

    # Verify they produce the same count
    if len(direct) != len(filtered):
        print(f"  ERROR: Count mismatch!")
    elif len(direct) != expected:
        print(f"  ERROR: Expected {expected}, got {len(direct)}")
    else:
        print(f"  ✓ Verified: {len(direct)} == {expected}")

print("\n" + "="*60)
print("CONCLUSION: Direct generation is much faster!")
print("="*60)
