from itertools import combinations
from sage.all import *


def _standardize(seq):
    order = {val: rank for rank, val in enumerate(sorted(seq), start=1)}
    return tuple(order[val] for val in seq)


def _avoids_patterns(perm, patterns):
    perm_list = list(perm)
    n = len(perm_list)
    for pattern_len, pattern_std in patterns:
        for positions in combinations(range(n), pattern_len):
            subseq = [perm_list[i] for i in positions]
            if _standardize(subseq) == pattern_std:
                return False
    return True

# Generate all Grassmannian permutations in S_n
# Method 1: Using pattern avoidance
def grassmannian_perms_patterns(n):
    """Generate Grassmannian perms by avoiding 321, 2143, 3142"""
    pattern_specs = tuple(
        (len(pattern), _standardize(pattern))
        for pattern in ((3, 2, 1), (2, 1, 4, 3), (3, 1, 4, 2))
    )
    return [p for p in Permutations(n) if _avoids_patterns(p, pattern_specs)]

# Method 2: Using descent condition (at most 1 descent)
def grassmannian_perms_descents(n):
    """Generate Grassmannian perms with at most one descent"""
    P = Permutations(n)
    return [p for p in P if len(p.descents()) <= 1]

# Example usage
n = 4
grass_perms = list(grassmannian_perms_patterns(n))
print(f"Grassmannian permutations in S_{n}: {len(grass_perms)} total")
print(f"Expected: {2**n - n} = {2**n - n}")

for p in grass_perms:
    print(f"{p} - descents at positions: {p.descents()}, inversions: {p.number_of_inversions()}")