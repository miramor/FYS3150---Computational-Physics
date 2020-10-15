import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import sys
import csv
import glob

#files = glob.glob('./Results/*.csv') #contains list of all files ending with csv, can be used if we wanna plot all files at once


systems = {
"systemA": ["Sun", "Earth"],
"systemB": ["Sun", "Earth", "Jupiter"],
"systemC": ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"],
"systemD": ["Sun", "Mercury", "Earth"]
}
MethodDic = {"E":"Euler", "VV": "Velocity Verlet", "EC": "Euler-Cromer"}
parameters = ['x_','y_','z_','vx_','vy_','vz_']

system = sys.argv[1]
method = sys.argv[2]
t_end = int(sys.argv[3])
h = float(sys.argv[4])
f_masses = open("masses.txt", 'r')
masses = [float(i) for i in f_masses.readlines()]
massesDic = {i:j for i,j in zip(systems["systemC"], masses)}



def plotEnergy(system):
    sys_names = [par + obj for obj in systems[system] for par in parameters]
    sys_data = pd.read_csv("Results/" + system + "_"+ method + ".csv", index_col=False, names=sys_names, skiprows=1)
    N = len(sys_data[sys_names[0]])





def plot_sys(system):
#     with open("test.csv", "r") as f:
#         reader = csv.reader(f)
#         row = next(reader)
#         print(row)

    sys_names = [par + obj for obj in systems[system] for par in parameters]
    sys_data = pd.read_csv("Results/" + system + "_"+ method + ".csv", index_col=False, names=sys_names, skiprows=1)
    N = len(sys_data[sys_names[0]])
    for i in range(int(len(sys_names)/6)):
        x = sys_data[sys_names[i*6]]
        y = sys_data[sys_names[i*6+1]]
        plt.plot(x,y, label = systems[system][i])

    """
    li = sys_data["x_Earth"].tail(1).index.item()
    xe0 = sys_data["x_Earth"][0]
    ye0 = sys_data["y_Earth"][0]
    ze0 = sys_data["z_Earth"][0]
    xel = sys_data["x_Earth"][li]
    yel = sys_data["y_Earth"][li]
    zel = sys_data["z_Earth"][li]

    absxe = abs(xe0-xel)
    absye = abs(ye0-yel)
    absze = abs(ze0-zel)
    #print(sys_data["x_Earth"][li])
    print(f"Check x position: {xe0}, {xel}, {absxe:.2e}")
    print(f"Check y position: {ye0}, {yel}, {absye:.2e}")
    print(f"Check z position: {ze0}, {zel}, {absze:.2e}")
    """

    plt.axis('equal')
    plt.legend()
    plt.ylabel("AU")
    plt.xlabel("AU")
    plt.title('%s - %s, h = %s'%(system, MethodDic[method],h))
    plt.savefig("Plots/" + system + "_" + method + "_" + f"{t_end}" + ".png", dpi=400)
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
    plt.savefig("Plots/" + system + "_" + method + "_" + f"{t_end}" + "_3D.png", dpi=400)
    plt.show()

plot_sys(system)
#plot3dPath(system)
