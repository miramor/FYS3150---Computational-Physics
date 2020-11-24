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

Two solve methods are created which call the Metropolis function for a set number of Monte Carlo cycles. The one called *solve_write()* writes continuous results for each sampling. This is useful for finding the effects of which way the lattice is initialized and how the spin-matrix and its observables evolve. The other one (*solve()*), just does not write the results to file for each sampling. However, a function called *writeFile()* allows to write the final results to file.

**How to use:**
Run *main.py* with the arguments in the commandline: N (lattice size) and MC_cycles (number of cycles). The user will be given an option to write the results for each sampling to file, looking only at one chosen temperature. The results are written to "Results/e_hist.csv" containing: mean energy per spin, mean absolute magnetization per spin, total number of accepted flips and energy per spin.
If option 2 is chosen (calculation for different temperatures in parallel), the user is then asked for start temperature, end temperature and number of datapoints wanted between the two temperatures (ideally 4, 8 or 12 etc to make use of max amount of threads). Then, it makes use of OpenMp and asks for number of threads you wish to use.  Finally the values T, \<E\>, \<M\>, Cv, chi are written to a .csv file named "Observable_[N]" in the "./Results" folder. The values and time needed for each temp is printed out in the terminal.
 
**Example use:**
```python
python3 main.py N MC_cycles
```

