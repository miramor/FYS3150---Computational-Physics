

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
n = 10
=======
print(pd.__version__)
n = 100
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb

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
<<<<<<< HEAD
#plt.savefig(f"Computed vs Exact sol n: {n}", dpi = 60)
plt.show()


#Make plot for time
=======

"""
#Add labels mid plot if wanted
n = x.size
plt.text(0.5, sol[int(n/2)] - 0.1, "Computed", fontsize=12, color="lightseagreen")
plt.text(0.5, exact[int(n/2)] + 0.1, "Exact", fontsize=12, color="crimson")
"""


#plt.savefig(f"Computed vs Exact sol n: {n}", dpi = 60)
plt.show()

>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb
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
<<<<<<< HEAD
plt.savefig("Timetosolve.png", dpi=60)
=======
#plt.savefig("Timetosolve.png", dpi=60)
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb

plt.show()




<<<<<<< HEAD
combinedTime = pd.concat([n, timeG, timeS], axis=1)
headers = ["n", "General time[s]", "Special time[s]"]
combinedTime.columns = headers
#vertical_stack.style.format("{:.2%}")
#print(vertical_stack)
vertical_stack.to_latex(index = False)
"""
=======

""" Maybe make a table is possible, seemed difficult. maybe use pandas plot
fig, ax = plt.subplots()
fig.set_figheight(15)
fig.set_figwidth(8)
ax.axis('off')
ax.axis('tight')
t= ax.table(cellText=dfTG.values, colWidths = [0.3]*len(dfTG.columns),  colLabels=dfTG.columns,  loc='center')
t.auto_set_font_size(False)
t.set_fontsize(12)
fig.tight_layout()
plt.show()
"""

fig,  ax = plt.subplots(1,1)
data1 = dfTG[["n", "time"]].head(6)
data2 = dfTS[["time"]].head(6)
vertical_stack = pd.concat([data1, data2], axis=1)
symbols = ["integration steps", "time general [sec]", "time special [sec]"]
vertical_stack.columns = symbols
#vertical_stack.style.format("{:.2%}")
#print(vertical_stack)

>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb
ax.axis("off")
pd.plotting.table(ax, vertical_stack, loc="upper right", colWidths=[0.3, 0.3, 0.3])
plt.savefig("Timetable.png")
plt.show()
<<<<<<< HEAD
"""
=======
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb
