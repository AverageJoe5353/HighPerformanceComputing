# Project 1: Optimization of the Matrix Product
Contributors: Addison, Connor, Joe, Rey

Report Authors: Addison, Joe, Rey

**revised-report.pdf** is the good copy of the report!

## Usage
1. Use `bash runBenchmarks.sh` to run the full set of tests for each language, compiler, optimization, and method.
2. Results will be saved in `results.csv`

Notes:
- The actual benchmarking script is in `benchmark.py` but I added the one-line shell script to meet the requirements
- Different values for N can be added to the `Ns` list in the configuration section of `benchmark.py`
- We used the newer intel compilers `icx` and `ifx`, if you use `icc` and `ifort` these will need to be replaced in the `compilers` dict in `benchmark.py`
