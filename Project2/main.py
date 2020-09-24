import os
import sys
import numpy as np


N = 3 #input("Choose n: ")
solution = "regular" #input("Choose solution: regular/potential: ")
"""
filename_plot = "_".join([solution, str(N)]) + ".pdf" #Name of figure file
filename_data = "_".join([solution, "N", str(N)]) + ".txt" #Name of data file.
plot_path = "/".join([".", "plots", solution]) #Directory to place the figure
data_path = "/".join([".", "results", solution]) #Directory to place the data file.
#First check if the directory exists. Otherwise, create it.
if not os.path.exists(data_path):
    os.makedirs(data_path) #Creates the directory
if not os.path.exists(plot_path):
    os.makedirs(plot_path)
"""

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output " + str(N) + " " + solution)

#os.system("echo creating plots...")
#os.system("python3 plot_results.py")

os.system("echo done.")

#Calculating analytical eigenvalues for V0 and V1
ana_eigvalV0 = np.zeros(N)
ana_eigvalV1 = np.zeros(N)
print("\n")
for i in range(0,N):
    #ana_eigvalV0[i] = (i+1)*np.pi**2
    #ana_eigvalV0[i] = (i+1)**2*np.pi**2
    ana_eigvalV1[i] = 3.0 + 4*i
#print("Analytic eigenvalues: V0\n", ana_eigvalV0)
print("Analytic eigenvalues: V1 \n", ana_eigvalV1)
