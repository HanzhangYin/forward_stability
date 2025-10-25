# Forward Stability Distribution

The original code is in `andy-s-original-code.sage`. The python version is`andy-s-original-code.py`.

The forward stability statistics are in `forward_stability_distribution.py`.
## Installation 
The code requires the following packages:
```bash
sage -pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org multipolynomial_bases tqdm
```

**Note:** The `--trusted-host` flags may be needed on some systems to bypass SSL certificate issues.

## Random Pair Selection Algorithm

### Two-Stage Sampling Strategy

### Stage 1: Pre-generate a Permutation Pool

**Purpose:** Create a manageable, diverse set of random permutations once.

**Implementation:**
```python
pool_size = min(POOL_SIZE, num_samples * 2)  # Default POOL_SIZE = 7000
perm_pool = [perms_obj.random_element() for _ in range(pool_size)]
```

**Key Points:**
- Pre-generates K random permutations (K = 7,000 by default)
- Each permutation is uniformly randomly selected from S_n using Sage's `random_element()`
- Pool is stored as a Python list to preserve Sage Permutation object methods
- Adaptive sizing: uses `min(POOL_SIZE, num_samples * 2)` to avoid generating more than needed

### Stage 2: Fast Index-Based Pairing

**Purpose:** Generate pairs efficiently using NumPy's vectorized operations.

**Implementation:**
```python
indices_u = np.random.randint(0, pool_size, num_samples)
indices_v = np.random.randint(0, pool_size, num_samples)
task_list = [(perm_pool[indices_u[i]], perm_pool[indices_v[i]]) 
             for i in range(num_samples)]
```

**Key Points:**
- Uses NumPy's `np.random.randint()` for batch index generation (~50-100× faster than loops)
- Generates num_samples pairs of indices in milliseconds
- Each index pair (i, j) is uniformly random from [0, pool_size)
- Pairs are sampled **with replacement** (same permutation can appear multiple times)

### Statistical Properties

**Sampling Space:**
- With pool_size = K, we sample uniformly from K² possible pairs
- Example: K = 7,000 → sampling from 49,000,000 possible pairs
- This is a large enough space to capture the statistical distribution

**Uniformity:**
- Each of the K² possible pairs has equal probability of selection: 1/K²
- The K permutations themselves are uniformly sampled from S_n
- Result: Approximately uniform sampling over a large representative subset

**Repetitions:**
- Pairs can repeat (sampling with replacement)
- Same permutation can be paired with itself (u = v is valid)
- Self-pairs where index i = j are allowed

### Performance Optimizations

1. **Pool Pre-generation:**
   - Generates K permutations once: O(K × n)
   - Avoids repeated Permutations(n) object creation

2. **NumPy Vectorization:**
   - Batch generates all indices: ~50-100× faster than Python loops
   - Pre-allocates arrays for memory efficiency

3. **List Storage:**
   - Keeps permutations as Sage objects (not NumPy arrays)
   - Preserves access to Sage methods like `.to_lehmer_code()`

4. **Multiprocessing:**
   - Distributes pair processing across all CPU cores
   - Uses `Pool.map()` with optimal chunking

### Configuration Parameters

```python
N = 110                    # S_N
NUM_RANDOM_SAMPLES = 100000  # Number of pairs to sample
POOL_SIZE = 7000           # Size of pre-generated permutation pool
```

**Tuning Guide:**
- **POOL_SIZE:** Larger = more diversity, slower setup. Recommended: 5,000-10,000
- **NUM_RANDOM_SAMPLES:** More samples = better statistical accuracy
- Optimal: `POOL_SIZE ≈ sqrt(NUM_RANDOM_SAMPLES)` to balance diversity vs. speed

### Algorithm Complexity

**Time Complexity:**
- Pool generation: O(K × n) where K = pool_size
- Index generation: O(num_samples) - vectorized, very fast
- Pair processing: O(num_samples × n²) - dominant cost, parallelized

**Space Complexity:**
- O(K × n) for storing the permutation pool
- O(num_samples) for storing indices and results


### Example Output

For n=110 with 100,000 samples:
```
Pre-generating pool of 7000 permutations...
Pool generation complete.
Generating 100000 random index pairs...
Building task list from random pairs...
Task list of 100000 pairs built.

Starting parallel processing with 12 cores...
Parallel processing finished in XXX seconds.
```

### Code Location

The random sampling implementation is in `forward_stability_distribution.py`, lines 191-205.

## Dependencies

- SageMath (for permutation operations)
- NumPy (for fast numerical operations)
- Python multiprocessing (for parallel execution)
- matplotlib (for visualization)
- tqdm (for progress bar display)
- multipolynomial_bases (for Schubert basis computations)
