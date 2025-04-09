"""
Riley Pantalena
CEE3610 - Hardy Cross
Hardy Cross Method for a Two-Loop Network with a Shared Branch
"""

import numpy as np

# Constants
N_EXPONENT = 1.852  # Exponent in the Hazen-Williams equation
MAX_ITERATIONS = 100  # Maximum number of iterations
TOLERANCE = 0.001  # Convergence tolerance

def kfunc(length, roughness, diameter):
    """
    Calculate the Hazen-Williams k constant given pipe attributes.

    Parameters:
        length (float): Pipe length in feet.
        roughness (float): Roughness coefficient (C value).
        diameter (float): Pipe diameter in inches.

    Returns:
        float: k such that head loss = k * |Q|^n * sign(Q)
    """
    return (4.727 * length) / (roughness ** 1.852 * diameter ** 4.8704)

def main():
    # Pipe parameters
    pipes = {
        1: {"length": 1500, "diameter": 10, "roughness": 100},
        2: {"length": 400, "diameter": 10, "roughness": 100},
        3: {"length": 600, "diameter": 10, "roughness": 120},
        4: {"length": 800, "diameter": 8, "roughness": 110},
        5: {"length": 700, "diameter": 8, "roughness": 110},
        6: {"length": 1000, "diameter": 10, "roughness": 90},
        7: {"length": 750, "diameter": 10, "roughness": 100}
    }
    
    # Calculate k values for each pipe
    k_values = {}
    for pipe_id, params in pipes.items():
        k_values[pipe_id] = kfunc(params["length"], params["roughness"], params["diameter"])
    
    # Define loops and initial flows
    # Loop 1: pipes 1, 2, 3, 4 (clockwise)
    # Loop 2: pipes 4, 5, 6, 7 (clockwise)
    
    # Initial flow guesses in cfs
    # Positive flow is in the direction specified by the loop orientation
    flows = {
        1: 0.77574,   # From A to D
        2: 2.0,       # From A to B
        3: 0.784386,  # From B to C
        4: 0.215614,  # From D to C (shared pipe)
        5: 0.963214,  # From E to D
        6: 1.431227,  # From E to F
        7: 0.215614   # From F to C
    }
    
    # Loop definitions: [pipe_id, direction]
    # Direction: 1 if flow is in same direction as loop traversal, -1 if opposite
    loop1 = [
        (1, 1),   # From A to D, positive in loop 1
        (4, 1),   # From D to C, positive in loop 1
        (3, -1),  # From C to B, negative in loop 1
        (2, -1)   # From B to A, negative in loop 1
    ]
    
    loop2 = [
        (5, 1),   # From E to D, positive in loop 2
        (4, -1),  # From C to D, negative in loop 2 (opposite of loop 1)
        (7, 1),   # From F to C, positive in loop 2
        (6, -1)   # From F to E, negative in loop 2
    ]
    
    # Perform Hardy-Cross iterations
    for iteration in range(MAX_ITERATIONS):
        # Calculate flow corrections for loop 1
        loop1_sum_hf = 0
        loop1_sum_dhf = 0
        
        for pipe_id, direction in loop1:
            flow = flows[pipe_id]
            k = k_values[pipe_id]
            
            # Head loss with direction: k * |Q|^n * sign(Q) * direction
            hf = k * abs(flow) ** N_EXPONENT * np.sign(flow) * direction
            loop1_sum_hf += hf
            
            # Derivative term: n * k * |Q|^(n-1)
            dhf = N_EXPONENT * k * abs(flow) ** (N_EXPONENT - 1)
            loop1_sum_dhf += dhf
        
        # Calculate correction for loop 1
        delta1 = -loop1_sum_hf / loop1_sum_dhf if loop1_sum_dhf != 0 else 0
        
        # Calculate flow corrections for loop 2
        loop2_sum_hf = 0
        loop2_sum_dhf = 0
        
        for pipe_id, direction in loop2:
            flow = flows[pipe_id]
            k = k_values[pipe_id]
            
            # Head loss with direction: k * |Q|^n * sign(Q) * direction
            hf = k * abs(flow) ** N_EXPONENT * np.sign(flow) * direction
            loop2_sum_hf += hf
            
            # Derivative term: n * k * |Q|^(n-1)
            dhf = N_EXPONENT * k * abs(flow) ** (N_EXPONENT - 1)
            loop2_sum_dhf += dhf
        
        # Calculate correction for loop 2
        delta2 = -loop2_sum_hf / loop2_sum_dhf if loop2_sum_dhf != 0 else 0
        
        # Apply corrections to the flows
        # For non-shared pipes in loop 1
        for pipe_id, direction in loop1:
            if pipe_id != 4:  # Non-shared pipes
                flows[pipe_id] += delta1 * direction
        
        # For non-shared pipes in loop 2
        for pipe_id, direction in loop2:
            if pipe_id != 4:  # Non-shared pipes
                flows[pipe_id] += delta2 * direction
        
        # For the shared pipe (pipe 4), we need to account for both loops
        # Loop 1 sees pipe 4 with direction 1, Loop 2 sees it with direction -1
        flows[4] += delta1 * 1 + delta2 * (-1)
        
        # Check for convergence
        if abs(delta1) < TOLERANCE and abs(delta2) < TOLERANCE:
            print(f"Converged after {iteration + 1} iterations")
            break
    else:
        print(f"Did not converge after {MAX_ITERATIONS} iterations")
    
    # Print final flows
    print("\nFinal Pipe Flows (cfs):")
    for pipe_id, flow in flows.items():
        print(f"Pipe {pipe_id}: {flow:.6f}")
    
    # Check mass balance at junctions
    junctions = {
        'A': {1: 1, 2: 1},         # Outflows from A
        'B': {2: -1, 3: 1},        # Inflow to B, outflow from B
        'C': {3: -1, 4: -1, 7: -1}, # Inflows to C
        'D': {1: -1, 4: 1, 5: -1},  # Inflow to D, outflow from D
        'E': {5: 1, 6: 1},         # Outflows from E
        'F': {6: -1, 7: 1}         # Inflow to F, outflow from F
    }
    
    print("\nJunction Mass Balance Check:")
    for junction, connections in junctions.items():
        balance = sum(flows[pipe_id] * direction for pipe_id, direction in connections.items())
        print(f"Junction {junction}: {balance:.6f}")
    
if __name__ == "__main__":
    main()
    
# =^.^=
#  |-|_/