"""
Riley Pantalena
CEE3610 - Hardy Cross
Coupled Hardy-Cross iteration for a two-loop network with a shared branch.
"""

import numpy as np


# Constants
N_EXPONENT = 1.852  # Exponent in the Hazen-Williams equation


def kfunc(length, roughness, diameter):
    """
    Calculate the Hazen-Williams k constant given pipe attributes.

    Parameters:
        length (float): Pipe length in feet.
        roughness (float): Roughness coefficient.
        diameter (float): Pipe diameter in inches.

    Returns:
        float: k such that head loss = k * |Q|^n (with sign attached separately).
    """
    return (4.727 * length) / (roughness ** 1.852 * diameter ** 4.8704)


def main():
    # Define loop 1: pipes 1, 2, 3, and shared pipe (pipe 4)
    # Pipe parameters: [length, diameter, roughness]
    # (pipes 1-3 are unique to loop 1, and pipe 4 is shared)
    k_values_1 = np.array([
        kfunc(1500, 100, 10),  # pipe 1
        kfunc(400, 100, 10),   # pipe 2
        kfunc(600, 120, 10),   # pipe 3
        kfunc(800, 110, 8)     # pipe 4 (shared)
    ])
    # Initial flow guesses in cfs for loop 1
    q1 = np.array([0.77574, 2.0, 0.784386, 0.2156135])
    # Assigned sign for head loss in loop 1 (defines pipe orientation)
    sign_1 = np.array([1, -1, -1, 1])

    # Define loop 2: pipes (shared pipe, then pipes 5, 6, 7)
    # For loop 2, the shared pipe appears with reversed orientation.
    k_values_2 = np.array([
        kfunc(800, 110, 8),   # pipe 4 (shared)
        kfunc(700, 110, 8),    # pipe 5
        kfunc(1000, 90, 10),   # pipe 6
        kfunc(750, 100, 10)    # pipe 7
    ])
    # Initial flow guesses in cfs for loop 2
    q2 = np.array([0.2156135, 0.9632135, 1.4312265, 0.2156135])
    # For loop 2, the sign on the shared branch is reversed compared to loop 1
    sign_2 = np.array([-1, -1, 1, 1])

    iterations = 100
    tolerance = 0.01

    # Iterative process:
    # For each loop, compute independent corrections for non-shared pipes, then update the shared pipe
    # by averaging the flows determined from each loop.
    for i in range(iterations):
        # ---- Loop 1: Compute head loss and derivative for each pipe ----
        h1 = np.zeros(4)
        dh1 = np.zeros(4)
        for j in range(4):
            # Head loss: k * |Q|^n, then apply the sign based on loop orientation
            h1[j] = k_values_1[j] * np.abs(q1[j]) ** N_EXPONENT * sign_1[j]
            # Derivative: n * k * |Q|^(n-1)
            dh1[j] = N_EXPONENT * k_values_1[j] * np.abs(q1[j]) ** (N_EXPONENT - 1)

        # In loop 1, pipes 1, 2, and 3 are non-shared (indices 0, 1, 2) and pipe 4 is shared.
        f1 = h1.sum()
        d1 = dh1.sum()
        delta1 = -f1 / d1  # Correction for loop 1

        # ---- Loop 2: Compute head loss and derivative for each pipe ----
        h2 = np.zeros(4)
        dh2 = np.zeros(4)
        for j in range(4):
            h2[j] = k_values_2[j] * np.abs(q2[j]) ** N_EXPONENT * sign_2[j]
            dh2[j] = N_EXPONENT * k_values_2[j] * np.abs(q2[j]) ** (N_EXPONENT - 1)

        # In loop 2, the shared pipe (index 0) is reversed.
        f2 = -h2[0] + h2[1] + h2[2] + h2[3]
        d2 = dh2.sum()
        delta2 = -f2 / d2  # Correction for loop 2

        # ---- Update flow rates (apply corrections) ----
        q1_new = q1.copy()
        q2_new = q2.copy()

        # For loop 1, update non-shared pipes (indices 0, 1, 2)
        q1_new[:3] = q1[:3] + delta1
        # For loop 2, update non-shared pipes (indices 1, 2, 3)
        q2_new[1:] = q2[1:] + delta2

        # ---- Update the shared pipe consistently ----
        # Loop 1 sees the shared pipe as q1[3] + delta1.
        # Loop 2 sees the shared pipe (reversed orientation) as q2[0] - delta2.
        # The new shared pipe flow is the average of the two.
        q_shared_new = ((q1[3] + delta1) + (q2[0] - delta2)) / 2
        q1_new[3] = q_shared_new
        q2_new[0] = q_shared_new

        # ---- Check for convergence ----
        if np.abs(delta1) < tolerance and np.abs(delta2) < tolerance:
            q1 = q1_new
            q2 = q2_new
            break

        q1 = q1_new
        q2 = q2_new

    # ----- Final Output -----
    print(f"Converged after {i+1} iterations")
    print("Final q values for loop 1 (in cfs):", q1)
    print("Final q values for loop 2 (in cfs):", q2)


if __name__ == "__main__":
    main()