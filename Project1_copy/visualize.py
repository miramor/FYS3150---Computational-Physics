

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
n = 10

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
#plt.savefig(f"Computed vs Exact sol n: {n}", dpi = 60)
plt.show()


#Make plot for time
dfTG = pd.read_csv("CPUtime general", names = ['n','time'])
dfTS= pd.read_csv("CPUtime special", names = ['n','time'])
print(dfTG["time"])
n = dfTG['n'].to_numpy()
timeG = dfTG['time'].to_numpy()
timeS = dfTS['time'].to_numpy()


plt.figure(figsize = (15, 8), dpi = 80)
plt.loglog(n, timeG, label = "General", color = "lightseagreen")
plt.loglog(n, timeS, linestyle = "dashed", label = "Special", color = "crimson")
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
plt.savefig("Timetosolve.png", dpi=60)

plt.show()




combinedTime = pd.concat([n, timeG, timeS], axis=1)
headers = ["n", "General time[s]", "Special time[s]"]
combinedTime.columns = headers
#vertical_stack.style.format("{:.2%}")
#print(vertical_stack)
vertical_stack.to_latex(index = False)
"""
ax.axis("off")
pd.plotting.table(ax, vertical_stack, loc="upper right", colWidths=[0.3, 0.3, 0.3])
plt.savefig("Timetable.png")
plt.show()
"""
