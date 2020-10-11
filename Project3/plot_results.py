import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd


def plot3dPath(x,y,z, navn):
    # Kan gjøre om til enten å få flere plott av ulike baner
    # eller flere i samme.
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z, label=f"{navn}")
    ax.legend()
    fig.savefig(f"Results/{navn}.png", dpi = 300)
    plt.show()

earth_data = pd.read_csv("Results/systemA.csv", index_col=False, names=['sx', 'sy', 'sz', 'svx', 'svy', 'svz','x', 'y', 'z', 'vx', 'vy', 'vz'])

x = earth_data['x'].to_numpy()
y = earth_data['y'].to_numpy()
z = earth_data['z'].to_numpy()

#print(earth_data.head())
plt.plot(x, y)
plt.axis('equal')
plt.savefig("sun_earth.png")


plot3dPath(x,y,z, "earth")
