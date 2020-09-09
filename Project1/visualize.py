import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

"""
How to use:

Simply run the visualize.py file and the showTables indicates
if tables are to be shown or not. Makes some plots for solutions
using Thomas algo vs exact.

This file was used to make tables using pandas and printing them.
Also to make it easy to copy straight into latex (using .to_latex)
"""

if not os.path.exists('ComputedvsExact'):
    os.makedirs('ComputedvsExact')

showTables = True #Prints all the tables if true.

n_list = [10, 10e1, 10e2, 10e3, 10e4, 10e5, 10e6] #used when making table of errors

for n in n_list[0:3]:
    dfG = pd.read_csv(f"./ResultsComputation/ResultsG_nval={int(n)}", names = ['x','solution', 'exact', 'log error'])
    dfS = pd.read_csv(f"./ResultsComputation/ResultsS_nval={int(n)}", names = ['x','solution', 'exact', 'log error'])

    x = dfG['x'].to_numpy()
    solG = dfG['solution'].to_numpy()
    solS = dfS['solution'].to_numpy()
    exact = dfG['exact'].to_numpy()

    plt.rcParams.update({'font.size' : 16})
    plt.figure(figsize = (15, 8), dpi = 80)
    plt.plot(x, solG, label = "Computed general", color = "lightseagreen")
    plt.plot(x, solS, linestyle = "dotted", label = "Computed special", color = "cyan", linewidth = 1)
    plt.plot(x, exact, linestyle = "dashed", label = "Exact", color = "crimson")
    plt.title(f"Computed vs Exact sol, n: {int(n)}", fontsize = 22)
    plt.xlabel("x_val")
    plt.ylabel("Solution (u)")
    plt.grid(axis='both', alpha=.22)

    plt.legend()

    # Remove borders
    plt.gca().spines["top"].set_alpha(0.0)
    plt.gca().spines["bottom"].set_alpha(0.6)
    plt.gca().spines["right"].set_alpha(0.0)
    plt.gca().spines["left"].set_alpha(0.6)
    plt.savefig(f"./ComputedvsExact/ComputedvsExact_sol_n{int(n)}.png", dpi = 60)
    #plt.show()

#--------------------------------------------------------------------------------------------------------

#Had to make new loop to avoid reading files where n is to big and causes memory error. Use only the relevant column and get max
error = []
for i in range(len(n_list)):
    dfG = pd.read_csv(f"./ResultsComputation/ResultsG_nval={int(n_list[i])}", names = ['x','solution', 'exact', 'log error'], usecols = ["log error"])
    error.append(dfG.max())
    #print(f"For n = {n_list[i]} Max error: {error[i]}")

#Make table with error values
dfLogErr = pd.DataFrame()
dfLogErr['n'] = n_list
dfLogErr['Max LogErr'] = error
#print(dfLogErr)
#print(dfLogErr.to_latex(index = False))

#-------------------------------------------------------------------------------------------------------------------------------------------------

#Table for 10 elements of the calculation results for the various methods.
dfG = pd.read_csv(f"./ResultsComputation/ResultsG_nval={10}", names = ['x','solution', 'exact', 'log error'])
dfS = pd.read_csv(f"./ResultsComputation/ResultsS_nval={10}", names = ['x','solution', 'exact', 'log error'])
LUsol = [0.4899 ,0.6119, 0.5986, 0.5356, 0.4542, 0.3660, 0.2754 ,0.1839, 0.0920] #Pulled from printing out solution in terminal from SolveLU function, missing deep precision
dfLU = pd.DataFrame()
dfLU['solution'] = LUsol
SolTable = pd.concat([dfG['x'], dfG['solution'], dfS['solution'], dfLU['solution'], dfG['exact']], axis=1, keys=['x_val', 'General', 'Special', 'LU', 'exact'])
#print(SolTable)
#print(SolTable.to_latex(index = False))


#--------------------------------------------------------------------------------------------------------

#Make plot for time and write out LATEX table.
#Read in all the time files and convert to numpy arrays to make it easy to plot
pd.set_option('display.float_format', '{:.3E}'.format) #Use scientific notation to display small numbers

dfTG = pd.read_csv("CPUtime general", names = ['n','time'])
dfTS = pd.read_csv("CPUtime special", names = ['n','time'])
dfLU = pd.read_csv("CPUtime LU", names = ['n','time'])

timeG = dfTG['time'].to_numpy()
n_G = dfTG['n'].to_numpy()
timeS = dfTS['time'].to_numpy()
n_S = dfTS['n'].to_numpy()
timeLU = dfLU['time'].to_numpy()
n_LU = dfLU['n'].to_numpy()

if timeG.size != timeS.size:#Just a check if calculated for all n with both models.
    print(f"Size of CPUtime general ({timeG.size}) and CPUtime special ({timeS.size}) does not match. Suggest running both for same n, if comparable results is wanted")



#Table of CPU times printed as latex table

dfTimeSpecial = dfTS['time']
dfTimeLU = dfLU['time']
#print(dfTimeSpecial)
combinedTime = pd.concat( [dfTG, dfTimeSpecial, dfTimeLU], axis = 1)
headers = ["n", "General time[s]", "Special time[s]", "LU time[s]"]
combinedTime.columns = headers
#print(combinedTime)
#print(combinedTime.to_latex(index = False))

if showTables:
    print("Table for error:")
    print(dfLogErr)
    print("\n \n")

    print("Table for solution for n = 10")
    print(SolTable)
    print("\n \n")

    print("Table for solution CPUtimes for different algorithms")
    print(combinedTime)




""" #Plotting of time for special and general, does not contain that interesting info
plt.figure(figsize = (15, 8), dpi = 80)
plt.plot(n_G, timeG, label = "General", color = "lightseagreen")
plt.plot(n_S, timeS, linestyle = "dashed", label = "Special", color = "crimson")
#plt.plot(n_LU, timeLU, linestyle = "dashed", label = "LUtime", color = "cyan") #Blows up the scale.
plt.title("Time to solve", fontsize = 16)
plt.grid(axis='both', alpha=.22)
plt.xlabel("Steps (n)")
plt.ylabel("Time [sec]")
plt.legend()

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.6)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.6)
plt.savefig(f"./Solvetime/Timetosolve.png", dpi=60)

#plt.show()
"""
