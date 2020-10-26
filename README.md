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

The two main parts of this program is a planet class and a solver class. In the planet class, we initialize planet objects. This object contains all the information about the planet, including name, mass, position and velocity. In the solver class, we have an Euler, Euler-Cromer and Velocity Verlet solver. The solver class also includes a function that writes position and velocity to a .csv file.

The solver function also has a test function that calculates the angular momentum and total energy in the first and the last step of the in the simulation, to see if these quantities are conserved.

We use *main.py* to compile all the .hpp and .cpp files into an executable, which then is executed automatically. *main.py* then runs *plot_results.py*.

**How to use:**
Run by writing "python3 main.py *system method end_time time_step*" in terminal.
The possible *systems* are systemA, systemB, systemC, and system E. Possible *methods* are E(Euler), EC(Euler-Cromer), VV(Velocity Verlet), VV2(Velocity Verlet used for perihelion precision). *end_time* is an integer choosing numbers of years to simulate. *time_step* is also given in the unit years. We used *time_step*=0.0001 for most of the simulations.

There will be instructions asking for inputs, these are always numbers(integers).

PS: We recommend plotting in 3D :) (axis != equal, which makes it look cooler)

**Run Example:**
main.py systemC VV 100 0.00001
*then use NASA data, adjust origin and plot 3D*
