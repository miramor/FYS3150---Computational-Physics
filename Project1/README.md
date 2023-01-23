### Project 1
Thomas algorithm and LU decomposition to solve the one-dimensional Poisson equation with Dirchlet boundary conditions. \
The Thomas alogrithm solves matrix equations of type Ax=b, where A is tridiagonal.
Hence the solution is computed with both a general algorithm and a specialised case, considering that the upper and lower diagonal are equal.
The algortihm includes forward and backward solution. The number of integration steps is denoted by *n*.

**How to use:**

Crreate a directory named **ResultsComputation**, which will contain all results of the run.
The C++ code can be run by using the included makefile
``` python
make all
```

Run make all by first setting *useSpecial = true* and then *useSpecial = false* in `main.cpp`. This does all the calculations and writes the results to the folder mentioned earlier. (To avoid excessive calculations, it is recommended to comment out the part which handles LU calculations, given in 2nd part of `main.cpp`, if it already was run before). Adjust the repeat paramater if needed (a low value will be compuationally faster with less precision, a high value will be slower and more precice). \
Dependencies: armadillo

When all results have been written for the various n values, type
``` python
python3 visualize.py
```
 into the terminal to display a table showing, among other things, errors and results for all methods for *n = 10*. In addition, some graphs are created and placed in a directoray to visualize the solution versus the exact curves.\
 Dependencies: os, numpy, matplotlib, pandas

