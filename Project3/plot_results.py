import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


earth_data = pd.read_csv("Results/systemA.csv", index_col=False, names=['sx', 'sy', 'sz', 'svx', 'svy', 'svz','x', 'y', 'z', 'vx', 'vy', 'vz'])

x = earth_data['x'].to_numpy()
y = earth_data['y'].to_numpy()

print(earth_data.head())
plt.plot(x, y)
plt.axis('equal')
plt.savefig("sun_earth.png")
