import sys
from collections import defaultdict

def get_step_probs(j):
    """
    Calculates the probabilities for the step Delta_j (where j = i).
    
    P(+1) = 1/j^2
    P(-1) = ((j-1)/j)^2
    P(0)  = 1 - P(+1) - P(-1)
    """
    if j < 2:
        # Delta_1 is not defined in this problem.
        return {0: 1.0}
        
    j = float(j)
    p_plus = 1.0 / (j * j)
    p_minus = ((j - 1.0) / j) ** 2
    
    # We calculate p_zero this way to avoid potential
    # floating-point precision errors from 2*(j-1)/(j*j)
    p_zero = 1.0 - p_plus - p_minus
    
    return {
        1: p_plus,
        0: p_zero,
        -1: p_minus
    }

def compute_E_M_i(i):
    """
    Computes E[M_i], which is the expected maximum of the random walk
    S_0, S_1, ..., S_{i-1}.
    
    This walk uses steps Delta_2, Delta_3, ..., Delta_i.
    
    Args:
        i (int): The (inclusive) index of the final step to use.
                 For E[M_3], i=3. For E[M_4], i=4.

    Returns:
        float: The expected maximum value.
    """
    
    # We start with the base case: M_1 = max(S_0) = max(0) = 0.
    # The distribution is: 100% probability of being at
    # position 0 with a max-so-far of 0.
    # Format: {(position, max_so_far): probability}
    dist = defaultdict(float)
    dist[(0, 0)] = 1.0

    # Loop through each step from Delta_2 up to Delta_i
    for j in range(2, i + 1):
        step_probs = get_step_probs(j)
        
        # Create a new distribution for the next state
        next_dist = defaultdict(float)
        
        # For every state in our current distribution...
        for (prev_pos, prev_max), prob in dist.items():
            
            # ...apply the 3 possible outcomes of the new step Delta_j
            for step_val, step_prob in step_probs.items():
                
                # 1. Calculate the new position
                new_pos = prev_pos + step_val
                
                # 2. Calculate the new maximum-so-far
                new_max = max(prev_max, new_pos)
                
                # 3. Add this path to the new distribution
                #    Its probability is the old path's prob * this step's prob
                next_dist[(new_pos, new_max)] += prob * step_prob
        
        # The new distribution becomes our current one for the next loop
        dist = next_dist

    # After the loop, 'dist' holds the final probability distribution
    # of all (position, max_so_far) pairs.
    
    # We get the final expectation by summing (max_so_far * probability)
    # over all states in our final distribution.
    total_expected_max = 0.0
    for (pos, max_val), prob in dist.items():
        total_expected_max += max_val * prob
        
    return total_expected_max

# --- Main execution ---
if __name__ == "__main__":
    
    # Set a higher recursion depth for the DP, just in case
    # (though this iterative DP doesn't use recursion)
    sys.setrecursionlimit(2000)

    print("--- Verifying the Agent's Calculations ---")
    
    # E[M_1] = E[max(S_0)] = E[max(0)] = 0
    E_M_1 = compute_E_M_i(1)
    print(f"E[M_1] (max after 0 steps): {E_M_1}")

    # E[M_2] = E[max(S_0, S_1)] = E[max(0, Delta_2)]
    # (1/4 * 1) + (1/2 * 0) + (1/4 * 0) = 1/4 = 0.25
    E_M_2 = compute_E_M_i(2)
    print(f"E[M_2] (max after 1 step):  {E_M_2}")

    # E[M_3] = E[max(S_0, S_1, S_2)]
    # This is the 1/3 calculation
    E_M_3 = compute_E_M_i(3)
    print(f"E[M_3] (max after 2 steps): {E_M_3:.10f} (Expected: 1/3 = {1/3:.10f})")

    # E[M_4] = E[max(S_0, S_1, S_2, S_3)]
    # This is the 13/36 calculation
    E_M_4 = compute_E_M_i(4)
    print(f"E[M_4] (max after 3 steps): {E_M_4:.10f} (Expected: 13/36 = {13/36:.10f})")

    print("\n--- Computing C_walk Lower Bound (by running for more steps) ---")
    # The lower bound E[M_i] will get closer and closer to the true value of C_walk
    # We can see how the value grows.
    
    E_M_10 = compute_E_M_i(10)
    print(f"E[M_10] (max after 9 steps):  {E_M_10:.10f}")
    
    E_M_20 = compute_E_M_i(20)
    print(f"E[M_20] (max after 19 steps): {E_M_20:.10f}")
    
    E_M_50 = compute_E_M_i(50)
    print(f"E[M_50] (max after 49 steps): {E_M_50:.10f}")
