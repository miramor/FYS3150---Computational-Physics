import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

systems = {
"systemA": ["Sun", "Earth"],
"systemB": ["Sun", "Earth", "Jupiter"],
"systemC": ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"],
"systemD": ["Sun", "Mercury", "Earth"]
}
parameters = ['x_','y_','z_','vx_','vy_','vz_']

def plot_sys(system):
    sys_names = [par + obj for obj in systems[system] for par in parameters]
    print(sys_names)
    sys_data = pd.read_csv("Results/" + system + ".csv", index_col=False, names=sys_names)
    for i in range(int(len(sys_names)/6)):
        x = sys_data[sys_names[i*6]]
        y = sys_data[sys_names[i*6+1]]
        plt.plot(x,y, label = systems[system][i])

    plt.axis('equal')
    plt.legend()
    plt.title('%s' %(system))
    plt.savefig("Plots/" + system + ".png", dpi=400)
    plt.show()


def plot3dPath(system):
    # Kan gjøre om til enten å få flere plott av ulike baner
    # eller flere i samme.
    sys_names = [par + obj for obj in systems[system] for par in parameters]
    print(sys_names)
    sys_data = pd.read_csv("Results/" + system + ".csv", index_col=False, names=sys_names)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i in range(int(len(sys_names)/6)):
        x = sys_data[sys_names[i*6]]
        y = sys_data[sys_names[i*6+1]]
        z = sys_data[sys_names[i*6+2]]
        ax.plot(x,y,z, label = systems[system][i])

    ax.legend()
    plt.savefig("Plots/" + system + "_3D.png", dpi=400)
    plt.axis('equal')
    plt.show()

plot3dPath("systemD")
