import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import sys
import csv
import math

systems = {
"systemA": ["Sun", "Earth"],
"systemB": ["Sun", "Earth", "Jupiter"],
"systemC": ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"],
"systemD": ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter"],
"systemE": ["Sun", "Mercury"]
}
MethodDic = {"E":"Euler", "VV": "Velocity Verlet", "EC": "Euler-Cromer", "VV2": "Velocity-Verlet"}
parameters = ['x_','y_','z_','vx_','vy_','vz_']

system = sys.argv[1]
method = sys.argv[2]
t_end = int(sys.argv[3])
h = float(sys.argv[4])
f_masses = open("masses.txt", 'r')
masses = [float(i) for i in f_masses.readlines()]
massesDic = {i:j for i,j in zip(systems["systemC"], masses)}


def plot_sys(system):
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
    plt.grid(axis='both', alpha=.22)
    # Remove borders
    plt.gca().spines["top"].set_alpha(0.0)
    plt.gca().spines["bottom"].set_alpha(0.6)
    plt.gca().spines["right"].set_alpha(0.0)
    plt.gca().spines["left"].set_alpha(0.6)
    plt.savefig("Plots/" + system + "_" + method + "_" + f"{t_end}" + ".pdf", dpi=200)
    plt.show()


def plot3dPath(system):
    # Kan gjøre om til enten å få flere plott av ulike baner
    # eller flere i samme.
    sys_names = [par + obj for obj in systems[system] for par in parameters]
    sys_data = pd.read_csv("Results/" + system + "_"+ method + ".csv", index_col=False, names=sys_names, skiprows=1)
    fig = plt.figure(figsize = (10,8))

    ax = fig.gca(projection='3d')
    N = len(sys_data[sys_names[0]])
    for i in range(int(len(sys_names)/6)):
        x = sys_data[sys_names[i*6]]
        y = sys_data[sys_names[i*6+1]]
        z = sys_data[sys_names[i*6+2]]
        ax.plot(x,y,z, label = systems[system][i])
    ax.legend(loc="center left", ncol=1)

    plt.title('%s - %s, h = %s' %(system, MethodDic[method],h))
    ax.set_zlabel("AU")
    ax.set_ylabel("AU")
    ax.set_xlabel("AU")
    plt.savefig("Plots/" + system + "_" + method + "_" + f"{t_end}" + "_3D.pdf", dpi=200)
    plt.show()

def plotOrbitDifference(filename, orbitform):
    """
    Plot distance to earth for system A and system B in order to see how Jupiter influences the Earth's orbit
    Must have list with csv files as input see definition of orbitA and orbitB below
    """
    sys_names_A = [par + obj for obj in systems["systemA"] for par in parameters]
    sys_names_B = [par + obj for obj in systems["systemB"] for par in parameters]
    if orbitform == "c":
        orbitA = pd.read_csv(filename[0], index_col=False, names=sys_names_A, skiprows=1)
        orbitB = pd.read_csv(filename[1], index_col=False, names=sys_names_B, skiprows=1)

    elif orbitform == "e":
        orbitA= pd.read_csv(filename[2], index_col=False, names=sys_names_A, skiprows=1)
        orbitB = pd.read_csv(filename[3], index_col=False, names=sys_names_B, skiprows=1)

    x_A = orbitA["x_Earth"]
    y_A = orbitA["y_Earth"]
    x_B = orbitB["x_Earth"]
    y_B = orbitB["y_Earth"]

    r_A = np.sqrt(x_A**2+y_A**2)
    r_B = np.sqrt(x_B**2+y_B**2)
    t_steg = np.linspace(0,12,120000) # h = 0.0001
    plt.plot(t_steg, r_B, label = "system B")
    plt.plot(t_steg, r_A, ":", label ="system A")
    plt.title("Earth's distance to Sun")
    plt.xlabel("time [year]")
    plt.ylabel("distance [AU]")
    plt.grid(axis='both', alpha=.22)
    plt.legend()
    # Remove borders
    plt.gca().spines["top"].set_alpha(0.0)
    plt.gca().spines["bottom"].set_alpha(0.6)
    plt.gca().spines["right"].set_alpha(0.0)
    plt.gca().spines["left"].set_alpha(0.6)

    plt.savefig("Plots/OrbitDifference" + "_" + f"{orbitform}" +".pdf", dpi=200)
    plt.show()

    plt.clf()
    plt.plot(t_steg, r_A-r_B)
    plt.title("Difference in Earth's distance to Sun")
    plt.grid()
    plt.ylim(-0.04,0.04)
    plt.xlabel('time [year]')
    plt.ylabel('distance [AU]')
    plt.savefig("Plots/OrbitDifferenceChange" + "_" + f"{orbitform}" +".pdf", dpi=200)


def calcAnglePerihelMerc():
    arcsecPerYr = 0.43 #arcseconds per year, 43'' per century
    #Convert to radians:
    radians = t_end*arcsecPerYr*math.pi/648000

    data = pd.read_csv(f"Results/Peri_Results.csv", index_col = False, names = ["x", "y", "z", "r", "i"])
    names = [par + obj for obj in systems["systemE"] for par in parameters]
    maxI =  data["x"].size-1 #max index
    MIO = data["i"].size-1 #index 0 - N, corresponds to last orbit perihelion, maxIndexOrbit

    # Position for perihelion in first and last orbit
    x_0 = data["x"][0]
    y_0 = data["y"][0]
    x_p = data["x"][maxI]
    y_p = data["y"][maxI]

    theta0 = math.atan(y_0/x_0)
    theta = math.atan(y_p/x_p)
    print(f"Angle: Numerical {theta:.4e} vs  Calculated {radians:.4e} after {t_end} years.  Theta0 {theta0:.2e}")
    print(f"Difference: {abs(theta-radians)}")

##filenames = ["Results/CE/systemA_VV_c12.csv", "Results/CE/systemB_VV_c12.csv", "Results/CE/systemA_VV_e12.csv","Results/CE/systemB_VV_e12.csv"]
#plotOrbitDifference(filenames, "c")
#plt.clf()
#plotOrbitDifference(filenames, "e")
#plt.clf()


if system == "systemE":
    choice = input("\nIf used NASA data press 0, otherwise the calculated theta is printed out:\n")
    if choice != "0":
        print(f"Calculates perihelion for {t_end} with h: {h}")
        calcAnglePerihelMerc()

else:
    #Always makes 2D plot:
    plot_sys(system)
    plt.clf()

    choice = input("Do you wanna plot in 3D?: (1)yes, (2)no:\n")
    if choice == "1":
        plot3dPath(system)
