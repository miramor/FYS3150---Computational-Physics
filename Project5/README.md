### Project 5

**Numerical Studies of the SIRS model (disease)**

This project makes use of both the fourth-order runge-kutta method (RK4) and monte-carlo sampling (MC) to simulate how a disease would spread among an isolated population consisting of three groups (states): susceptible (S), infected (I) and recovered (R). We study differnt rates of recovery that decide how (fast) the sick population recovers (I->R). In addition, the basic SIRS model (std) is extended by adding vital dynamics (VD), which includes deaths and births. Moreover, vaccines (Vac) allow a direct a direct transition from the susceptible state to the recovered one, reducing the impact of the diesease. Furthermore, an option exist to add seasonal variation (SV) to the rate of transmission (a), which affects the transition from susceptible to sick. The rate of immunity loss (c) is set constant.

**Structure:**
The program is divided into one big class named `SIRS.cpp` (+SIRS.hpp) which contains all the solver methods for the RK4 and Monte-carlo method and writing to file. Each instance of the class can solve for the given timesteps for either MC or RK4. The choice of solver type and problem type is done with distinct specialization methods which are called upon after the constructor (such as specRK4_VD). `main.cpp` includes all the paramaters needed for all the different runs. It creates an object for both RK4 and MC, and, depending on choices made in the terminal by the user, it will make solve the given problem for the four different recovery-values (b) and writes them to their designated folder in the folder named *Results* as .csv files. 

The .csv files contains one header line including information of the problem solved which is used exctract data about the specific run when plotting. The columns after are as follows: S, I, R. Once its all finished writing to file we use the function *Plot_HealthsStatus* in `plot.py` to visualize the values. is plotted and placed in the correct folders(being PlotsVac, Plotsstd and PlotsVD). 
The `main.py` files initially runs the makefile which compiles all the files. Then executes the code and then runs `plot.py` finalizing the program.

**Extra:**
Code not used when running includes a `plot2.py`, which looked to increase options but currently gives some errors. When it comes to the deadthInfected.pdf plot it was a modified version of this, but to avoid major bugs or changes to a functioning program a 2nd folder was made; named Project5_copy. This program loops over mulitple dI(deaths disease) or f(vaccine) and writes the final death toll to a file which then was plotted using the *plot_deathInfected* method. 

**How to use:**
Simply run the python program `main.py` and make decisions on if to enable certain options such as: seasonal variation, vital dynamics and vaccines by typing 0 (no) or 1 (yes).
```python
python3 main.py
```

**Example use:**
Using seasonal variation combined with vital dynamics, plotted relevant figures but no relative difference & expecation values printed out. The text in the paranthesis explains what each choice is:
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
// 0, 1, * -> using vital dynamics, no seasonal
// 1, 0, 0 -> using standard, with seasonal
```
**Dependencies** \
gcc: armadillo, lapack, blas \
python3: seaborn, matplotlib, pandas, numpy

