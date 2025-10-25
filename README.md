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

- Uses NumPy's `np.random.randint()` for batch index generation (~50-100× faster than loops)
- Each index pair (i, j) is uniformly random from [0, pool_size)
- Pairs are sampled **with replacement** (same permutation can appear multiple times)


### List Storage:
   - Keeps permutations as Sage objects (not NumPy arrays)
   - Preserves access to Sage methods like `.to_lehmer_code()`

### Multiprocessing:
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

## Generating the Chart

To generate the forward stability distribution chart, run:

```bash
sage forward_stability_distribution.py
```

**What it does:**
1. Samples `NUM_RANDOM_SAMPLES` random pairs (u, v) from S_n × S_n
2. Computes FS(u,v) = max(integer_support(u,v)) + 1 for each pair
4. Creates a bar chart saved as `sN_stability_sampled_SAMPLES.png`

**Output includes:**
- Expected value E[FS(u,v)]
- Standard deviation (σ) and variance (σ²)

**Progress tracking:** The script displays a real-time progress bar showing completion percentage, processing speed, and estimated time remaining.

## Dependencies

- SageMath (for permutation operations)
- NumPy (for fast numerical operations)
- Python multiprocessing (for parallel execution)
- matplotlib (for visualization)
- tqdm (for progress bar display)
- multipolynomial_bases (for Schubert basis computations)
