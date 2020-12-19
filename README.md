# FYS3150
FYS3150 – Computational Physics .  Filer til undervisning og oppgaver. Prosjekter.

## Projects
### Project 1:
Thomas algorithm and LU decomposition to solve second order differential equation. Writing the results to files for different n values and visualizing results.

**How to use:**

Create a directory called **ResultsComputation** in project 1 folder, this houses all the results.
Run make all for useSpecial = true once and then again with  useSpecial = false. This does all the calculations and writes the results to the folder mentioned earlier. (Can comment out the part which handles LU calculations, 2nd part of main, if already ran main.cpp once to avoid excessive calculations). Adjust repeat paramater if needed (low = fast and less precision *** high = slow and more precice time). 

When all results for the different n values are written. Then run visualize.py if you wanna see tables which shows things such as error and results for all methods for *n = 10*. Also creates some plots which is placed in a directoray to visualize the solution vs exact curves.

### Project 2:
Using the Jacobi algorithm to solve eigenvalue problems.

**How to use:**
We defined the potentials V0, V1 and V2, respectively for the buckling beam, quantum dots for one electron, and quantum dots for two electrons. Run main.py and follow instructions. This scripts runs the makefile which compiles all the .cpp and .hpp files into one executable file. After that is done, it automatically runs the executable file which then solves the problems and and writes the results to text files. Thereafter, the plot_results.py script reads the data and save plots as .pdf's. All text files are moved to the "results" folder, and all plots are moved to the "plots" folder.

In order to use the function that records the CPU run time and number of iterations for the Buckling Beam problem, uncomment the function "writeTimeResults(omega)" in the main.cpp file, and the two functions at the bottom of "plot_results.py".

In order to run the tests, uncomment the test functions in main.cpp.


### Project 3:
**Simulating the solar system**

The two main parts of this program is a planet class and a solver class. In the planet class, we initialize planet objects. This object contains all the information about the planet, including name, mass, position and velocity. In the solver class, we have an Euler, Euler-Cromer and Velocity Verlet solver. The solver class also includes a function that writes position and velocity to a .csv file. In addition for the perihelion part (perihelion precession of Mercury), theres a special method which only writes the perihelion data which uses modified version of the standard methods used for the other solvers. For this and the VV method, the estimated time for the whole calculation and progress will be printed out.

The solver function also has a test function that calculates the relative angular momentum and total energy in the first and the last step of the in the simulation, to see if these quantities are conserved.

We use *main.py* to compile all the .hpp and .cpp files into an executable, which then is executed automatically. *main.py* then runs *plot_results.py* which currently plots in 2d, 3d and calculates the angle used when looking at the perihelion shift.

**How to use:**
Run by writing "python3 main.py *system method end_time time_step*" in terminal.
The possible *systems* are systemA (Sun-Earth), systemB (Sun-Earth-Jupiter), systemC (Sun plus planets of solar system), and systemE (Sun-Mercury). Possible *methods* are E (Euler), EC (Euler-Cromer), VV (Velocity Verlet), VV2 (Velocity Verlet only used for perihelion preseccion of Mercury). *end_time* is an integer choosing numbers of years to simulate. *time_step* is also given in the unit years. We used *time_step*=0.0001 for most of the simulations, while for the perihelion for 1 year 1e-7 or 1e-8 seemed sufficient. 

There will be instructions asking for inputs, these are always numbers(integers). 
These question will be about: asking to plot3d, adjusting center of mass, using custom/NASA data when reading in initial data for the planets. Some combinations of methods and solvers are not compatible (such as systemA and VV2), for which the program will auto-fix itself and give the user feedback of what it has changed.

Make sure the message "*methodName* started" (VV started) is printed before leaving the terminal alone, to avoid wasting time.

PS: We recommend plotting in 3D :) (axis != equal, which makes it look cooler)

**Run Examples:**

1.Plots all planets excluding Pluto in 3D. Calculations takes about 2 min.
main.py systemC VV 100 0.00001
*then use NASA data, adjust origin and plot 3D*

2.For perihelion:
This *h* = 1e-8 takes about 50s.
main.py systemE VV2 1 0.00000001
*will then automatically choose the right method*
*last input, dont press 0 when asked: If used NASA data press 0, otherwise the calculated theta is printed out*

3.Basic example
main.py systemA VV 2 0.00001

### Project 4:
**Studies of phase transitions in magnetic systems**
The isingModel make use of the Metropolis algorithm (Markov chain Monte Carlo) to find the most stable state for the system depending on the size of the lattice (N) and the temperature (T) given. The program picks a random place in the grid, calculates the change in energy a spin-filp would cause, and depending on the transition probability, the spin is either flipped or remains the same. This sampling is repeated for N^2 x MCcycles. Therby, system specific properties such as energy, magnetization, specific heat capacity and susceptibility are computed.

The program is quite simply built up of an isingModel class and a *main.cpp* which set ups the parallelization over different temperatures using OpenMP. Additionally a program *plot.py* is used to visualize the results. The *main.py* simply compiles and executes all the code. The solve method loops over N^2 x MCcycles of iterations and calls on the Metropolis function which either flips the spin of a given position in the lattice or keeps its orientation. 

Two solve methods are created which call the Metropolis function for a set number of Monte Carlo cycles. The one called *solve_write()* writes continuous results for each sampling. This is useful for finding the effects of which way the lattice is initialized and how the spin-matrix and its observables evolve. The other one (*solve()*), just does not write the results to file for each sampling. However, a function called *writeFile()* allows to write the final results to file. Note that there is a variable cutoff, which stands for the number of Monte Carlo cycles after which the mean values are collected. 

**How to use:**
Run *main.py* with the arguments in the commandline: N (lattice size) and MC_cycles (number of cycles). The user will be given an option (*Option 1*) to write the results for each sampling to file, looking only at one chosen temperature. The results are written to "Results/e_hist.csv" containing: mean energy per spin, mean absolute magnetization per spin, total number of accepted flips and energy per spin.
If *Option 2* is chosen (calculation for different temperatures in parallel for a given lattice), the user is then asked for start temperature, end temperature and number of datapoints wanted between the two temperatures (ideally 4, 8 or 16 etc to make use of max amount of threads depending on CPU). Then, it makes use of OpenMp and asks for number of threads you wish to use, pick a number between 1 and the max num of threads which is given in the terminal (if invalid pick, auto changes to max threads).  Finally the values T, \<E\>, \<|M|\>, Cv, chi are written to a .csv file named "Observable_[N]" in the "./Results" folder. The values and time needed for each temp is printed out in the terminal.

Using *Option 2* a folder named **"Results"** is added for the programs which takes hours to run. In our case, the observable files were computed using 15 million cycles for N={40,60,80,100}. Currently *plot.py* uses these files. If plots for stabilization development are needed run *Option 1* with e.g 400 cycles and N=20. Then run *plot_stabi()* function at line 247 in *plot.py*. Similarly, a histogram for the probability distribution can be plotted by uncommenting the *plot_hist_nump()* function. We now used 100 000 cycles instead. (You may need to change T in title, annotation location depends on which temp used and savefig name.)

 
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

### Project 5:
**Numerical Studies of the SIRS model (disease)**
This project makes use of both the fourth-order runge-kutta method (RK4) and monte-carlo sampling (MC) to simulate how a disease would spread among an isolated population consisting of three groups (states): susceptible (S), infected (I) and recovered (R). We study differnt rates of recovery that decide how (fast) the sick population recovers (I->R). In addition, the basic SIRS model (std) is extended by adding vital dynamics (VD), which includes deaths and births. Moreover, vaccines (Vac) allow a direct a direct transition from the susceptible state to the recovered one, reducing the impact of the diesease. Furthermore, an option exist to add seasonal variation (SV) to the rate of transmission (a), which affects the transition from susceptible to sick. The rate of immunity loss (c) is set constant.

**Structure:**
The program is divided into one big class named SIRS.cpp (+SIRS.hpp) which contains all the solver methods for the RK4 and Monte-carlo method and writing to file. Each instance of the class can solve for the given timesteps for either MC or RK4. The choice of solver type and problem type is done with distinct specialization methods which are called upon after the constructor (such as specRK4_VD). *main.cpp* includes all the paramaters needed for all the different runs. It creates an object for both RK4 and MC, and, depending on choices made in the terminal by the user, it will make solve the given problem for the four different recovery-values (b) and writes them to their designated folder in the folder named *Results* as .csv files. 

The .csv files contains one header line including information of the problem solved which is used exctract data about the specific run when plotting. The columns after are as follows: S, I, R. Once its all finished writing to file we use the function *Plot_HealthsStatus* in *plot.py* to visualize the values. is plotted and placed in the correct folders(being PlotsVac, Plotsstd and PlotsVD). 
The *main.py* files initially runs the makefile which compiles all the files. Then executes the code and then runs *plot.py* finalizing the program.

**Extra:**
Code not used when running includes a plot2.py, which looked to increase options but currently gives some errors. When it comes to the deadthInfected.pdf plot it was a modified version of this, but to avoid major bugs or changes to a functioning program a 2nd folder was made; named Project5_copy. This program loops over mulitple dI(deaths disease) or f(vaccine) and writes the final death toll to a file which then was plotted using the *plot_deathInfected* method. 

**How to use:**
Simply run the python program main.py and make decisions on if to enable certain options such as: seasonal variation, vital dynamics and vaccines by typing 0(no) or 1(yes).
```python
python3 main.py
```

**Example use:**
Using seasonal variation combined with vital dynamics, plotted relevant figures but no relative difference & expecation values printed out. Used () to explain what each choice is:
```python
python3 main.py
1 (use seasonal variation)
1 (use vital dynamics)
(main.cpp runs, plot.py gets run)
0 (do not print out extra data like error and expecation values)
```

There's a maximum of 4 in total, but just 3 inputs needed when for the the main.cpp file (last one is for plot.py). Each combination resulting in a different type of problem. Some example combinations is noted below (* = meaning no input here):
```python
// 0, 0, 0 -> using standard
// 1, 0, 1 -> using vaccines with seasonal variation
// 1, 1, * -> using vital dynamics, with seasonal
// 0, 0, 1 -> using vaccines without seasonal variation
// 0, 1, * -> vital dynamics, no seasonal
// 1, 0, 0 -> using standard, with seasonal
```


