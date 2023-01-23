### Project 4:
Studies of phase transitions in magnetic systems

The isingModel makes use of the Metropolis algorithm (Markov chain Monte Carlo) to find the most stable state for the system depending on the size of the lattice (N) and the temperature (T) given. The program picks a random place in the grid, calculates the change in energy a spin-filp would cause, and depending on the transition probability, the spin is either flipped or remains the same. This sampling is repeated for N^2 x MCcycles. Therby, system specific properties such as energy, magnetization, specific heat capacity and susceptibility are computed.

The program is quite simply built up of an isingModel class and a *main.cpp* which set ups the parallelization over different temperatures using OpenMP. Additionally a program *plot.py* is used to visualize the results. The *main.py* simply compiles and executes all the code. The solve method loops over N^2 x MCcycles of iterations and calls on the Metropolis function which either flips the spin of a given position in the lattice or keeps its orientation. 

Two solve methods are created which call the Metropolis function for a set number of Monte Carlo cycles. The one called *solve_write()* writes continuous results for each sampling. This is useful for finding the effects of which way the lattice is initialized and how the spin-matrix and its observables evolve. The other one (*solve()*), just does not write the results to file for each sampling. However, a function called *writeFile()* allows to write the final results to file. Note that there is a variable cutoff, which stands for the number of Monte Carlo cycles after which the mean values are collected. 

**How to use:**
Run `main.py` with the arguments in the commandline: N (lattice size) and MC_cycles (number of cycles). The user will be given an option (*Option 1*) to write the results for each sampling to file, looking only at one chosen temperature. The results are written to "Results/e_hist.csv" containing: mean energy per spin, mean absolute magnetization per spin, total number of accepted flips and energy per spin.
If *Option 2* is chosen (calculation for different temperatures in parallel for a given lattice), the user is then asked for start temperature, end temperature and number of datapoints wanted between the two temperatures (ideally 4, 8 or 16 etc to make use of max amount of threads depending on CPU). Then, it makes use of OpenMp and asks for number of threads you wish to use, pick a number between 1 and the max num of threads which is given in the terminal (if invalid pick, auto changes to max threads).  Finally the values T, \<E\>, \<|M|\>, Cv, chi are written to a .csv file named "Observable_[N]" in the "./Results" folder. The values and time needed for each temp is printed out in the terminal.

Using *Option 2* a folder named **"Results"** is added for the programs which takes hours to run. In our case, the observable files were computed using 15 million cycles for N={40,60,80,100}. Currently *plot.py* uses these files. If plots for stabilization development are needed run *Option 1* with e.g 400 cycles and N=20. Then run *plot_stabi()* function at line 247 in `plot.py`. Similarly, a histogram for the probability distribution can be plotted by uncommenting the *plot_hist_nump()* function. We now used 100 000 cycles instead. (You may need to change T in title, annotation location depends on which temp used and savefig name.)

 
**Example use:**
Option 1a: Run time with i7-4790k with N=20, 100 000 cycles and T=2.4 (using one thread). Use *plot_hist_nump()* to visualize this result.
```python
python3 main.py 20 100000
1
2.4

```

Option 1b: Run time with i7-4790k with N=20, 400 cycles and T=2.4 (using one thread). Use *plot_stabi()* to visualize this result.
```python
python3 main.py 20 400
1
2.4
```

Option 2: Computed files to be used by the functions "T_critical()" and "plot_all_obs(\[40,60,80,100\])" in *plot.py*. Run time with i7-4790k with 5 million cycles on 8 threads is 790-850 seconds.
```python
python3 main.py 40 5000000
2
1 2.4 8
8
```
**Dependencies** \
gcc: - \
python3: seaborn, matplotlib, pandas, numpy, scipy, sklearn