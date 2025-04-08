"""
Riley Pantalena
CEE3610 - Hardy Cross
"""

import numpy as np

# Constants
n = 1.852  # Exponent in Hazen-Williams equation

# Hazen-Williams k function: length in feet, diameter in inches, Q in cfs, headloss in feet
def kfunc(L, C, D):
    return ((4.727 * L)) / (C ** (1.852) * D ** (4.8704))

# Length (L) in feet, Diameter (D) in inches, and flow rates (Q) in cfs (no conversion to MGD)
kvalues1 = np.array([kfunc(1500, 100, 10),  # L = 1500 feet, D = 10 inches
                     kfunc(400, 100, 10),
                     kfunc(600, 120, 10),
                     kfunc(800, 110, 8)])  # D = 8 inches

Qvalues1 = np.array([0.77574,  # Flow rates in cfs (no conversion to MGD)
                     2.0,
                     0.784386,
                     0.2156135])

sign1 = np.array([1, -1, -1, 1])

kvalues2 = np.array([kfunc(800, 110, 8), 
                     kfunc(700, 110, 8),
                     kfunc(1000, 90, 10),
                     kfunc(750, 100, 10)])

Qvalues2 = np.array([0.2156135,
                     0.9632135,
                     1.4312265,
                     0.2156135])

sign2 = np.array([-1, -1, 1, 1])

iterations = 100
tolerance = 0.01
iteration = 0

# Iterative solution loop
for i in range(iterations):
    # Loop 1 headloss and flow calculation
    headloss1 = kvalues1 * np.abs(Qvalues1) ** n  # Headloss in feet
    headlosssign1 = headloss1 * sign1

    # Calculate delta Q for loop 1
    hloverQ1 = headlosssign1 / Qvalues1
    deltaQ1 = (-headlosssign1.sum()) / (n * np.abs(hloverQ1.sum()))

    # Update Qvalues1
    Qvalues1_new = Qvalues1 + deltaQ1

    # Loop 2 headloss and flow calculation
    headloss2 = kvalues2 * np.abs(Qvalues2) ** n  # Headloss in feet
    headlosssign2 = headloss2 * sign2

    # Calculate delta Q for loop 2
    hloverQ2 = headlosssign2 / Qvalues2
    deltaQ2 = (-headlosssign2.sum()) / (n * np.abs(hloverQ2.sum()))

    # Update Qvalues2
    Qvalues2_new = Qvalues2 + deltaQ2

    # Average Q for pipe 4 (shared between loops)
    # Compute headloss for the shared pipe as contributed by both loops
    H_shared = (kvalues1[3] * abs(Qvalues1[3])**n) * sign1[3] + (kvalues2[0] * abs(Qvalues2[0])**n) * sign2[0]

    # Compute the derivative contributions (dH/dQ) for the shared pipe from both loops
    dH_dQ_shared = n * ( (kvalues1[3] * abs(Qvalues1[3])**(n-1)) + (kvalues2[0] * abs(Qvalues2[0])**(n-1)) )

    # Then update the shared pipe with a single correction
    deltaQ_shared = - H_shared / dH_dQ_shared
    Q_shared_new = (Qvalues1[3] + deltaQ_shared + Qvalues2[0] + deltaQ_shared) / 2

    # Set this corrected flow for both loop representations
    Qvalues1_new[3] = Q_shared_new
    Qvalues2_new[0] = Q_shared_new

    # Check convergence
    if abs(deltaQ1) < tolerance and abs(deltaQ2) < tolerance:
        break
    Qvalues1 = Qvalues1_new
    Qvalues2 = Qvalues2_new

    iteration += 1

# Output the results
print(f"Converged after {iteration} iterations")
print(deltaQ1)
print(deltaQ2)
print("Final Qvalues1 (in cfs):", Qvalues1)
print("Final Qvalues2 (in cfs):", Qvalues2)
