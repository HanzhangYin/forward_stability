#!/usr/bin/env sage -python
"""
Test script to demonstrate the progress bar during Grassmannian permutation generation.

This script will:
1. Clear any cached results for the test value of n
2. Generate Grassmannian permutations from scratch (showing progress bar)
3. Run a second time to demonstrate cache loading (no progress bar)
"""

import sys
import os
import shutil

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sage.all import *
import time

# Import the functions from our module
from grassmannian_stability_distribution import (
    get_grassmannian_perms,
    CACHE_DIR
)

# Test with n=15 (32,753 permutations - enough to see progress bar)
TEST_N = 15

print("="*70)
print("GRASSMANNIAN PERMUTATION GENERATION - PROGRESS BAR TEST")
print("="*70)
print()

# Clear cache for this n value
cache_file = os.path.join(CACHE_DIR, f"grass_n{TEST_N}.pkl")
if os.path.exists(cache_file):
    print(f"Clearing cache for n={TEST_N} to demonstrate generation...")
    os.remove(cache_file)
    print()

# Test 1: Generate from scratch (should show progress bar)
print(f"TEST 1: First run (generate from scratch)")
print("-" * 70)
start = time.time()
perms1 = get_grassmannian_perms(TEST_N, use_cache=True)
elapsed1 = time.time() - start
print(f"Total time: {elapsed1:.3f} seconds")
print()

# Test 2: Load from cache (should be instant, no progress bar)
print(f"TEST 2: Second run (load from cache)")
print("-" * 70)
start = time.time()
perms2 = get_grassmannian_perms(TEST_N, use_cache=True)
elapsed2 = time.time() - start
print(f"Total time: {elapsed2:.3f} seconds")
print()

# Summary
print("="*70)
print("SUMMARY")
print("="*70)
print(f"Permutations generated: {len(perms1):,}")
print(f"Expected (2^{TEST_N} - {TEST_N}): {2**TEST_N - TEST_N:,}")
print(f"First run (with generation): {elapsed1:.3f}s")
print(f"Second run (from cache): {elapsed2:.3f}s")
print(f"Speedup: {elapsed1/elapsed2:.1f}x faster")
print()
print("If you saw a progress bar during TEST 1, everything is working!")
print("="*70)
