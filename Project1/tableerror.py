import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

n_vals = [10, 10e1, 10e2] #add more values when calculated
dict = {}
errorlist = []
nlist = []

#fetch first errorvalue from each file, each column
for n in n_vals:
    # Read in the first line of each results for the different n values and just use the log error column
    df = pd.read_csv(f"Results nval = {int(n)}", names = ['x','solution', 'exact', 'log error'], usecols = ["log error"], nrows = 1)
    error = df["log error"].iloc[0]
    errorlist.append(error)

dfError = pd.DataFrame({"n":n_vals,"log error":errorlist})

dfError.to_csv("test.csv")

#Make table for logerror
fig,  ax = plt.subplots(1,1,figsize=(5,2))
data1 = dfError[["n", "log error"]]
#vertical_stack.style.format("{:.2%}")
ax.axis("off")
pd.plotting.table(ax, data1, loc="upper right", colWidths=[0.3, 0.3, 0.3])
plt.savefig("LogError.png")
plt.show()
