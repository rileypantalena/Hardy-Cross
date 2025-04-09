# Hardy-Cross Method for Pipe Networks

## Overview

This Python implementation solves pipe network flow distribution problems using the Hardy-Cross method. It specifically handles two-loop networks with a shared pipe, calculating flow rates in each pipe while maintaining both head loss balance around loops and mass conservation at junctions.

## Background

The Hardy-Cross method is an iterative technique for analyzing flow distribution in closed-loop pipe networks. It was developed by Hardy Cross in 1936 and remains a fundamental approach in hydraulic engineering for solving water distribution networks.

The method works by:
1. Making initial flow estimates that satisfy continuity at each junction
2. Calculating head losses around each loop
3. Computing flow corrections to balance energy (head loss) around each loop
4. Iteratively applying these corrections until convergence

## Features

- Calculates pipe resistance coefficients using the Hazen-Williams formula
- Handles multi-loop networks with shared pipes
- Verifies mass balance at network junctions
- Provides detailed output of flow distribution
- Monitors convergence with user-defined tolerance

## Usage

### Input Parameters

Define your pipe network with the following information:

```python
# Pipe parameters
pipes = {
    pipe_id: {
        "length": length_in_feet,
        "diameter": diameter_in_inches,
        "roughness": hazen_williams_c_value
    }
}

# Define loops as lists of (pipe_id, direction) tuples
# Direction: 1 if flow direction matches loop traversal, -1 if opposite
loop1 = [(pipe_id, direction), ...]
loop2 = [(pipe_id, direction), ...]

# Initial flow estimates (must satisfy continuity at junctions)
flows = {pipe_id: initial_flow_estimate, ...}
```

### Example Network

This implementation was created for a network with the following structure:

```
A -> 1 -> D <- 5 <- E
|         |         |
V         V         V
2         4         6
|         |         |
V         V         V
B -> 3 -> C <- 7 <- F
```

Where:
- Loop 1 consists of pipes 1, 4, 3, and 2
- Loop 2 consists of pipes 5, 4, 7, and 6
- Pipe 4 is shared between both loops

### Running the Code

```bash
python hardy_cross.py
```

## Configuration Parameters

- `N_EXPONENT`: The exponent used in the Hazen-Williams formula (typically 1.852)
- `MAX_ITERATIONS`: Maximum number of iterations before terminating if convergence is not achieved
- `TOLERANCE`: Convergence criterion for flow corrections

## Output

The program outputs:
1. Number of iterations required for convergence
2. Final flow rates in each pipe
3. Mass balance verification at each junction

## Mathematical Background

The Hardy-Cross method uses the relationship between flow rate and head loss:

`h = k * |Q|^n * sign(Q)`

Where:
- `h` is the head loss
- `k` is the resistance coefficient
- `Q` is the flow rate
- `n` is the flow exponent (1.852 for Hazen-Williams)

Flow corrections are calculated as:

`ΔQ = -Σ(h) / Σ(n * k * |Q|^(n-1))`

## Handling Shared Pipes

For shared pipes, corrections from each loop are applied directly according to the pipe's orientation in each loop:

`Q_shared += ΔQ_loop1 * direction_in_loop1 + ΔQ_loop2 * direction_in_loop2`

This approach maintains consistency with both the conservation of energy principle and continuity requirements.