import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

n_vals = [10, 10e1, 10e2] #add more values when calculated
dict = {}
dfN = pd.DataFrame(n_vals)
dfError = pd.DataFrame([])

for n in n_vals:
    # Read in the first line of each results for the different n values and just use the log error column
    df = pd.read_csv(f"Results nval = {int(n)}", names = ['x','solution', 'exact', 'log error'], usecols = ["log error"], nrows = 1)
    error = [df["log error"].values[0]]
    dfError.append(error)

print(dfError)


#Make table for logerror
fig,  ax = plt.subplots(1,1)
data1 = df[["x", "log error"]].head(6)
#vertical_stack.style.format("{:.2%}")
ax.axis("off")
pd.plotting.table(ax, data1, loc="upper right", colWidths=[0.3, 0.3, 0.3])
plt.savefig("LogError.png")
plt.show()
