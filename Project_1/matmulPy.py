import time
import random
import numpy as np

def threeloop(N):
    print(f"Running Python three-loop matrix multiplication for N={N}...")
    rows = N
    cols = N    
    A = [[random.uniform(0, 10) for i in range(cols)] for j in range(rows)] # Fill matrix A with random floats between 0 and 10
    B = [[random.uniform(0, 10) for i in range(cols)] for j in range(rows)] # Fill matrix B with random floats between 0 and 10
    C = [[0.0 for i in range(cols)] for j in range(rows)] # Initialize result matrix to zeros

    start = time.time()
    for i in range(rows):
        for j in range(cols):
            C[i][j] = 0.0       # Set entry ij to zero (not strictly necessary here)
            for k in range(rows): # Dot Prod of A_i and B_j
                C[i][j] += A[i][k] * B[k][j]
    end = time.time()
    wall_time = end - start
    return float(wall_time)

def npmatmul(N):
    print(f"Running Python numpy matrix multiplication for N={N}...")
    rows = N
    cols = N
    A = np.random.uniform(0, 10, (rows, cols)) # Fill matrix A with random floats between 0 and 10
    B = np.random.uniform(0, 10, (rows, cols)) # Fill matrix B with random floats between 0 and 10

    start = time.time()
    C = np.matmul(A, B)
    end = time.time()
    wall_time = end - start
    return float(wall_time)