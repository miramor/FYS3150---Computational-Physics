### Project 2:
Implementation of the Jacobi alogrithm to solve second-order eigenvalue equations.

**How to use:**
We defined the potentials V0, V1 and V2, respectively for the buckling beam, quantum dots for one electron, and quantum dots for two electrons. 
Run the code by typing
``` python
 python3 main.py
```
and follow instructions. This scripts runs the makefile which compiles all the `.cpp` and `.hpp` files into one executable file. After that is done, it automatically runs the executable file which then solves the problems and and writes the results to text files. Thereafter, the `plot_results.py` script reads the data and save plots as .pdf's. All text files are moved to the "results" folder, and all plots are moved to the "plots" folder.

In order to use the function that records the CPU run time and number of iterations for the Buckling Beam problem, uncomment the function "writeTimeResults(omega)" in the `main.cpp` file, and the two functions at the bottom of `plot_results.py`.

In order to run the tests, uncomment the test functions in `main.cpp`.

**Dependencies** \
gcc: armadillo, lapack, blas \
python3: numpy, matplotlib, sys, pandas