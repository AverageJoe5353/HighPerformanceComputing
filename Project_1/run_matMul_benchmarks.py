import os
import random
import time
import numpy as np


################################ Configuration ################################
Ns = [512, 1024, 2048, 4096] # Matrix sizes to test

langs = ['C', 'Fortran', 'Python'] # List languages here
compilers = {
    'C': ['gcc','icx'], # List C compilers here
    'Fortran': ['gfortran','ifx'], # List Fortran compilers here
    'Python': ['Python'] # Python does not require a compiler
}
methods = {
    'C': ['Three-Loop'], # Methods used in C
    'Fortran': ['Three-Loop'], # Methods used in Fortran
    'Python': ['Three-Loop', 'NumPy'] # Methods used in Python
}
###############################################################################


#### Functions to run tests for each language and method ####

# C Three-Loop Method
def runCtest(n, optimization, compiler="gcc"):
    exe_name = f"matmulC"
    compile_cmd = f"{compiler} -o {exe_name} matmulC.c -O{optimization} -lm"
    os.system(compile_cmd)
    print(f"Running C three-loop matrix multiplication for N={n} with -O{optimization}...")
    run_cmd = f"./{exe_name} {n}"
    wall_time = os.popen(run_cmd).read().strip()
    return float(wall_time)

# Fortran Three-Loop Method
def runFortrantest(n, optimization, compiler="gfortran"):
    exe_name = f"matmulFortran"
    compile_cmd = f"{compiler} -o {exe_name} matmulFortran.f90 -O{optimization} -lm"
    os.system(compile_cmd)
    print(f"Running Fortran three-loop matrix multiplication for N={n} with -O{optimization}...")
    run_cmd = f"./{exe_name} {n}"
    wall_time = os.popen(run_cmd).read().strip()
    return float(wall_time)

# Python Three-Loop Method
def runNaivePytest(N):
    print(f"Running Python three-loop matrix multiplication for N={N}...")
    rows = N
    cols = N    
    A = [[random.uniform(0, 10) for i in range(cols)] for j in range(rows)]
    B = [[random.uniform(0, 10) for i in range(cols)] for j in range(rows)]
    C = [[0.0 for i in range(cols)] for j in range(rows)]

    start = time.time()
    for i in range(rows):
        for j in range(cols):
            for k in range(rows):
                C[i][j] += A[i][k] * B[k][j]
    end = time.time()
    wall_time = end - start
    return float(wall_time)

# Python NumPy Method
def runNumPytest(N):
    print(f"Running Python numpy matrix multiplication for N={N}...")
    rows = N
    cols = N
    A = np.random.uniform(0, 10, (rows, cols))
    B = np.random.uniform(0, 10, (rows, cols))

    start = time.time()
    C = np.matmul(A, B)
    end = time.time()
    wall_time = end - start
    return float(wall_time)


#### Loop over languages, compilers, optimization levels, and methods to run tests ####
CPU_model = os.popen("lscpu | grep 'Model name' | sed 's/Model name:\\s*//'").read().strip()
with open("results.csv", "w") as f:
    f.write(f"OS,CPU,Language,Compiler,-O Level,Method,N={Ns[0]},N={Ns[1]},N={Ns[2]},N={Ns[3]}\n")
    
    for lang in langs:
        for method in methods[lang]:
            # Python Test (No Compiler or Optimization)
            if lang == 'Python':
                wall_times = []
                for N in Ns:
                    if method == 'Three-Loop':
                        wall_time = runNaivePytest(N)
                    elif method == 'NumPy':
                        wall_time = runNumPytest(N)
                    wall_times.append(wall_time)
                f.write(f"Linux,{CPU_model},{lang},Python,-O0,{method}," + ",".join(f"{wt}" for wt in wall_times) + "\n")
            else:
                for compiler in compilers[lang]:
                    for optimization in [0, 1, 2, 3]:
                        wall_times = []
                        for N in Ns:
                            if lang == 'C':
                                wall_time = runCtest(N, optimization, compiler)
                            elif lang == 'Fortran':
                                wall_time = runFortrantest(N, optimization, compiler)
                            wall_times.append(wall_time)
                        f.write(f"Linux,{CPU_model},{lang},{compiler},-O{optimization},{method}," + ",".join(f"{wt}" for wt in wall_times) + "\n")
    


