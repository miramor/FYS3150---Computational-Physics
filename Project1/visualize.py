

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#n_choice = input("Choose a n values to visualize: (10, 10e1, 10e2, 10e3)\n")
#n = int(n_choice)
#n_list = [10, 10e1, 10e2, 10e3, 10e4, 10e5]
n_list = [10, 10e1]
for n in n_list:
    dfG = pd.read_csv(f"./ResultsComputation/ResultsG_nval={int(n)}", names = ['x','solution', 'exact', 'log error'])
    dfS = pd.read_csv(f"./ResultsComputation/ResultsS_nval={int(n)}", names = ['x','solution', 'exact', 'log error'])

    x = dfG['x'].to_numpy()
    solG = dfG['solution'].to_numpy()
    solS = dfG['solution'].to_numpy()
    exact = dfG['exact'].to_numpy()


    plt.figure(figsize = (15, 8), dpi = 80)
    plt.plot(x, solG, label = "Computed general", color = "lightseagreen")
    plt.plot(x, solS, linestyle = "dotted", label = "Computed special", color = "cyan", linewidth = 0.4)
    plt.plot(x, exact, linestyle = "dashed", label = "Exact", color = "crimson")
    plt.title(f"Computed vs Exact sol, n: {n}", fontsize = 16)
    plt.xlabel("steps (n)")
    plt.ylabel("Solution (u)")
    plt.grid(axis='both', alpha=.22)
    plt.legend()

    # Remove borders
    plt.gca().spines["top"].set_alpha(0.0)
    plt.gca().spines["bottom"].set_alpha(0.6)
    plt.gca().spines["right"].set_alpha(0.0)
    plt.gca().spines["left"].set_alpha(0.6)
    plt.savefig(f"./ComputedvsExact/ComputedvsExact_sol_n{n}.png", dpi = 60)
    plt.show()


    #Make plot for time
    dfTG = pd.read_csv("CPUtime general", names = ['n','time'])
    dfTS= pd.read_csv("CPUtime special", names = ['n','time'])
    dfLU= pd.read_csv("CPUtime LU", names = ['n','time'])
    #print(dfTG["time"])
    n = dfTG['n'].to_numpy()
    timeG = dfTG['time'].to_numpy()
    timeS = dfTS['time'].to_numpy()
    timeLU = dfLU['time'].to_numpy()

    if timeS.size != timeG.size:
        print("Has to be same nr of values in CPUtime general and CPUtime special! Exiting program")
        break


    plt.figure(figsize = (15, 8), dpi = 80)
    plt.plot(n, timeG, label = "General", color = "lightseagreen")
    plt.plot(n, timeS, linestyle = "dashed", label = "Special", color = "crimson")
    #plt.plot(n, timeLU, linestyle = "dashed", label = "LUtime", color = "cyan")
    plt.title("Time to solve", fontsize = 16)
    plt.grid(axis='both', alpha=.22)
    plt.xlabel("steps (n)")
    plt.ylabel("time [sec]")
    plt.legend()

    # Remove borders
    plt.gca().spines["top"].set_alpha(0.0)
    plt.gca().spines["bottom"].set_alpha(0.6)
    plt.gca().spines["right"].set_alpha(0.0)
    plt.gca().spines["left"].set_alpha(0.6)
    plt.savefig(f"./Solvetime/Timetosolve.png", dpi=60)

    plt.show()



    #Include this if want a Latex formatted table of time assuming the files
    # CPUtime* files exists.

    combinedTime = pd.concat([n, timeG, timeS, timeLU], axis=1)
    headers = ["n", "General time[s]", "Special time[s], "LU time[s]""]
    combinedTime.columns = headers
    vertical_stack.to_latex(index = False)
    
