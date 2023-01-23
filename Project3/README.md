### Project 3:
Simulating the solar system

The two main parts of this program are a planet class and a solver class. In the planet class, we initialize planet objects. This object contains all the information about the planet, including name, mass, position and velocity. In the solver class, we have an Euler, Euler-Cromer and Velocity Verlet solver. The solver class also includes a function that writes position and velocity to a .csv file. In addition for the perihelion part (perihelion precession of Mercury), theres a special method which only writes the perihelion data which uses modified version of the standard methods used for the other solvers. For this and the VV method, the estimated time for the whole calculation and progress will be printed out.

The solver function also has a test function that calculates the relative angular momentum and total energy in the first and the last step of the in the simulation, to see if these quantities are conserved.

We use `main.py` to compile all the `.hpp` and `.cpp` files into an executable, which then is executed automatically. `main.py` then runs `plot_results.py` which currently plots in 2d, 3d and calculates the angle used when looking at the perihelion shift.

**How to use:**
Run by writing 
``` python
 python3 main.py <system method> <end_time time_step>
```
into the terminal.
The possible *systems* are systemA (Sun-Earth), systemB (Sun-Earth-Jupiter), systemC (Sun plus planets of solar system), and systemE (Sun-Mercury). Possible *methods* are E (Euler), EC (Euler-Cromer), VV (Velocity Verlet), VV2 (Velocity Verlet only used for perihelion preseccion of Mercury). *end_time* is an integer choosing numbers of years to simulate. *time_step* is also given in the unit years. We used *time_step*=0.0001 for most of the simulations, while for the perihelion for 1 year 1e-7 or 1e-8 seemed sufficient. 

There will be instructions asking for inputs, these are always numbers (integers). 
These question will be about: asking to plot3d, adjusting center of mass, using custom/NASA data when reading in initial data for the planets. Some combinations of methods and solvers are not compatible (such as systemA and VV2), for which the program will auto-fix itself and give the user feedback of what it has changed.

Make sure the message "*methodName* started" (VV started) is printed before leaving the terminal alone, to avoid wasting time.

PS: We recommend plotting in 3D :) (axis != equal, which makes it look cooler)

**Run Examples:**

1. Plots all planets excluding Pluto in 3D. Calculations takes about two minutes.
``` python
 python3 main.py systemC VV 100 0.00001
```
then use NASA data, adjust origin and plot 3D.

2. For perihelion:
Setting *h* = 1e-8 the code runs for about 50s.
``` python
python3 main.py systemE VV2 1 0.00000001
```
The correct method will automatically be chosen.
For the last input, do not enter 0 when asked: If used NASA data press 0, otherwise the calculated theta is printed out*

3. Basic example
``` python
python3 main.py systemA VV 2 0.00001
```

**Dependencies** \
gcc: armadillo \
python3: numpy, matplotlib, sys, pandas, math