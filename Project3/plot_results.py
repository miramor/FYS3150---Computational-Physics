import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import sys
import csv
import glob

files = glob.glob('./Results/*.csv')
print(files)

systems = {
"systemA": ["Sun", "Earth"],
"systemB": ["Sun", "Earth", "Jupiter"],
"systemC": ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"],
"systemD": ["Sun", "Mercury", "Earth"]
}
MethodDic = {"E":"Euler", "VV": "Velocity Verlet", "EC": "Euler-Cromer"}
parameters = ['x_','y_','z_','vx_','vy_','vz_']

system = sys.argv[1]
print(system)

h = int(sys.argv[4])

method = sys.argv[3]
print(method)

# f = open("Results/" + system + ".csv", "r")
# method = f.readline().strip()
# f.close()


def plot_sys(system):
    with open("test.csv", "r") as f:
        reader = csv.reader(f)
        row = next(reader)
        print(row)

    sys_names = [par + obj for obj in systems[system] for par in parameters]
    sys_data = pd.read_csv("Results/" + system + "_"+ method + ".csv", index_col=False, names=sys_names, skiprows=1)
    N = len(sys_data[sys_names[0]])
    for i in range(int(len(sys_names)/6)):
        x = sys_data[sys_names[i*6]]
        y = sys_data[sys_names[i*6+1]]
        plt.plot(x,y, label = systems[system][i])

    plt.axis('equal')
    plt.legend()
    plt.ylabel("AU")
    plt.xlabel("AU")
    plt.title('%s - %s, h = %s'%(system, MethodDic[method],h))

    plt.savefig("Plots/" + system + "_" + method + "_" + h + ".png", dpi=400)
    plt.show()



def plot3dPath(system):
    # Kan gjøre om til enten å få flere plott av ulike baner
    # eller flere i samme.
    sys_names = [par + obj for obj in systems[system] for par in parameters]
    sys_data = pd.read_csv("Results/" + system + "_"+ method + ".csv", index_col=False, names=sys_names, skiprows=1)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    N = len(sys_data[sys_names[0]])
    for i in range(int(len(sys_names)/6)):
        x = sys_data[sys_names[i*6]]
        y = sys_data[sys_names[i*6+1]]
        z = sys_data[sys_names[i*6+2]]
        ax.plot(x,y,z, label = systems[system][i])

    ax.legend()
    plt.title('%s - %s, h = %s' %(system, MethodDic[method],h))
    ax.set_zlabel("AU")
    ax.set_ylabel("AU")
    ax.set_xlabel("AU")
    plt.savefig("Plots/" + system + "_" + method + "_" + h + "_3D.png", dpi=400)
    plt.show()
plot_sys(system)
#plot3dPath(system)
