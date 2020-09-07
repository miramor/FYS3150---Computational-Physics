

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

n_choice = input("Choose a n values to visualize: (10,100, 1 000, 10 000, 100 000)\n")
n = int(n_choice)

df = pd.read_csv(f"Results nval = {n}", names = ['x','solution', 'exact', 'log error'])

x = df['x'].to_numpy()
sol = df['solution'].to_numpy()
exact = df['exact'].to_numpy()


plt.figure(figsize = (15, 8), dpi = 80)
plt.plot(x, sol, label = "Computed", color = "lightseagreen")
plt.plot(x, exact, linestyle = "dashed", label = "Exact", color = "crimson")
plt.title(f"Computed vs Exact sol, n: {n}", fontsize = 16)
plt.grid(axis='both', alpha=.22)
plt.legend()

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.6)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.6)
plt.savefig(f"./ComputedvsExact/Computed vs Exact sol n: {n}", dpi = 60)
plt.show()


#Make plot for time
dfTG = pd.read_csv("CPUtime general", names = ['n','time', 'timeArma'])
dfTS= pd.read_csv("CPUtime special", names = ['n','time', 'timeArma'])
#print(dfTG["time"])
n_val = dfTG['n'].to_numpy()
timeG = dfTG['time'].to_numpy()
timeS = dfTS['time'].to_numpy()
timeArma_spec = dfTG['timeArma'].to_numpy()
timeArma_gen = dfTS['timeArma'].to_numpy()


plt.figure(figsize = (15, 8), dpi = 80)
plt.plot(n_val, timeG, label = "General", color = "lightseagreen")
plt.plot(n_val, timeS, label = "Special", color = "crimson")
#plt.plot(n_val, timeArma_spec, linestyle = "dashed",label = "General Arma", color = "lightseagreen")
#plt.plot(n_val, timeArma_spec, linestyle = "dashed", label = "Special Arma", color = "crimson")
plt.title("Time to solve", fontsize = 16)
plt.grid(axis='both', alpha=.22)
plt.xlabel("integration steps")
plt.ylabel("time [sec]")
plt.legend()

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.6)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.6)
plt.savefig(f"./Solvetime/Timetosolve_n_{n}.png", dpi=60)

plt.show()



#Include this if want a Latex formatted table of time assuming the files
# CPUtime* files exists.

combinedTime = pd.concat([n, timeG, timeS], axis=1)
headers = ["n", "General time[s]", "Special time[s]"]
combinedTime.columns = headers
vertical_stack.to_latex(index = False)
