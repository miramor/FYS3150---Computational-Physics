import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import sys
import csv
import math

#files = glob.glob('./Results/*.csv') #contains list of all files ending with csv, can be used if we wanna plot all files at once


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
    print(sys_names)
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

#plot_sys(system)
#plt.clf()

#plot3dPath(system)
#plt.clf()

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

##filenames = ["Results/CE/systemA_VV_c12.csv", "Results/CE/systemB_VV_c12.csv", "Results/CE/systemA_VV_e12.csv","Results/CE/systemB_VV_e12.csv"]
#plotOrbitDifference(filenames, "c")
#plt.clf()
#plotOrbitDifference(filenames, "e")
#plt.clf()
def calcAnglePerihelMerc():
    arcsecPerYr = 0.43 #arcseconds per year, 43'' per century

    #Convert to radians:
    radians = arcsecPerYr*math.pi/648000

    data = pd.read_csv(f"Results/perihelioMerc.csv", index_col = False, names = ["x", "y", "r", "i"])

    names = [par + obj for obj in systems["systemE"] for par in parameters]
    #print(names)
    #planetData = pd.read_csv("Results/systemE_VV.csv", index_col = False, names = names,  usecols = ["x_Sun", "y_Sun"] ,skiprows=1) #, usecols = ["x_Sun", "y_Sun", "z_Sun"]


    maxI =  data["x"].size-1 #max index
    MIO = data["i"].size-1 #index 0 - N, corresponds to last orbit perihelion, maxIndexOrbit
    print("x0 xn, y0 yn")
    #print(planetData["x_Sun"][0])
    #print(planetData["y_Sun"][0])
    #print(planetData["x_Sun"][MIO])
    #print(planetData["y_Sun"][MIO])
    #print(planetData["z_Sun"][MIO])

    #sunDrift_x = planetData["x_Sun"][MIO]
    #sunDrift_y = planetData["y_Sun"][MIO]

    #print(f"Drift Sun: {sunDrift_x} ,  {sunDrift_y}")



    # Position for perihelion in first and last orbit
    x_0 = data["x"][0]
    y_0 = data["y"][0]
    x_p = data["x"][maxI]
    y_p = data["y"][maxI]

    #Correct for movement of sun
    #y_p += -sunDrift_x
    #x_p += -sunDrift_y
    theta0 = math.atan(y_0/x_0)
    theta = math.atan(y_p/x_p)
    print(f"y_p/x_p --->   {y_p/x_p}")
    print(f"Angle: Numerical {theta:.4e} vs  Calculated {radians} after {t_end} years.  Theta0 {theta0}")


plot_sys(system)
calcAnglePerihelMerc()

"""
1 0.0000001
Ingen kraft: Numerical -3.1531e-05 vs  Calculated 2.0846988287710047e-06 after 1 years.
Med kraft:   Numerical -3.1531e-05 vs  Calculated 2.0846988287710047e-06 after 1 years.

"""

"""
#Started making function to plot first and last orbit for mercury
def checkPerihelAngleMerc():
    sys_names = [par + obj for obj in systems["systemE"] for par in parameters]
    data = pd.read_csv(f"Results/systemE_VV.csv", index_col=False, names=sys_names, skiprows=1)
    #print(data["x_Sun"])
    N = int(t_end/h)
    print(N)
    ISLO = int(N - (95/(365*t_end))/h) #IndexStartLastOrbit, index to make sure we measure the right value
    #x_min = data["x_Mercury"][ISLO]
    #y_min = data["y_Mercury"][ISLO]
    print(int(ISLO))
    x_0 = data["x_Mercury"][0]
    y_0 = data["y_Mercury"][0]
    merc_x = data["x_Mercury"][ISLO:N]
    merc_y = data["y_Mercury"][ISLO:N]

    x_max = merc_x.max()
    y_max = merc_y.max()
    xi_max = merc_x.idxmin()
    yi_max = merc_y.idxmin()

    print(f"Max x = {x_max} for index {xi_max}")
    print(f"Max x = {y_max} for index {yi_max}")
"""
#checkPerihelAngleMerc()
