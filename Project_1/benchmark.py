import os
import matmulPy # Import the Python matrix multiplication module

################################ Configuration ################################
Ns = [512, 1024, 2048, 4096] # Matrix sizes to test

# Ns = [64, 128, 256, 512] # Reduced sizes for debugging

langs = ['C', 'Fortran', 'Python'] # List languages here
compilers = {
    'C': ['gcc'], # List C compilers here
    'Fortran': ['gfortran'], # List Fortran compilers here
    'Python': ['Python'] # Python does not require a compiler
}
methods = {
    'C': ['Three-Loop'], # Methods used in C
    'Fortran': ['Three-Loop'], # Methods used in Fortran
    'Python': ['Three-Loop', 'NumPy'] # Methods used in Python
}
###############################################################################


#### Function to compile and execute C and Fortran code ####
def time_threeloop(n, lang, optimization, compiler):
    exe_name = f"threeLoop_{lang}" # Name of the executable to be created
    ext = 'c' if lang == 'C' else 'f90' # File extension for C or Fortran
    compile_cmd = f"{compiler} -o {exe_name} matmul{lang}_3Loop.{ext} -O{optimization} -lm" # Compilation command
    os.system(compile_cmd) # Compile the C code
    print(f"Running {lang} three-loop matrix multiplication for N={n} with -O{optimization}...") # Print status
    run_cmd = f"./{exe_name} {n}" # Command to run the executable
    wall_time = os.popen(run_cmd).read().strip() # Execute and capture output
    return float(wall_time)



#### Loop over languages, compilers, optimization levels, and methods to run tests ####
with open("results.csv", "w") as f:
    f.write(f"OS,Language,Compiler,-O Level,Method," + ",".join(f"N={N}" for N in Ns) + "\n")
    
    for lang in langs:
        for method in methods[lang]:
            # Special case for Python since it does not require compilation or have optimization levels
            if lang == 'Python': 
                wall_times = []
                for N in Ns:
                    if method == 'Three-Loop':
                        wall_time = matmulPy.threeloop(N) # Executes the three-loop method
                    elif method == 'NumPy':
                        wall_time = matmulPy.npmatmul(N) # Executes the NumPy method
                    wall_times.append(wall_time)
                f.write(f"Linux,{lang},Python,-O0,{method}," + ",".join(f"{wt}" for wt in wall_times) + "\n")
            else:
                for compiler in compilers[lang]:
                    for optimization in [0, 1, 2, 3]:
                        wall_times = []
                        for N in Ns:
                            wall_time = time_threeloop(N, lang, optimization, compiler) # Run test using the generic function
                            wall_times.append(wall_time)
                        f.write(f"Linux,{lang},{compiler},-O{optimization},{method}," + ",".join(f"{wt}" for wt in wall_times) + "\n")
    


