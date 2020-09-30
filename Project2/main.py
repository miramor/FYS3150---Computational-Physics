import os
import sys
import numpy as np
"""
run main.py and follow instructions.
"""

n = int(input("Choose n value (10 - 200, delta n = 10):  "))
solution = input("Choose potential type (V0, V1 or V2):  ").upper()
omega = "0"
dictRho = {"V0" : "1", "V1" : "4.6", "V2" : "9"}
rho_max = dictRho[solution]

if solution == "V2":
    dict = {"1" : "0.25", "2" : str(1/20), "3" : str(1/54.7386)}
    omegaChoice = input("Choose omega value( 1 = 1/4, 2 = 1/20, 3 = 1/54.7386 ):  ")
    omega_rho = {"0.25": "10", str(1/20): "20", str(1/54.7386): "30"}
    omega = dict[omegaChoice]
    rho_max = omega_rho[omega]
    print(omega)

#Defining the command line argument used when executing main.cpp and plot_results.py
cla = str(n) + " " + solution + " " + rho_max + " " + omega


#Making paths for results and plots
if solution == "V2":
    filename_data = "_n_" + str(n) + "_" + solution + "_w_%.6f.txt" %(float(omega))
    filename_plot = "eigplot_" + solution + "_n_" + str(n) + "_w_%.3f.PDF" %(float(omega))
else:
    filename_data = "_n_" + str(n) + "_" + solution + ".txt"
    filename_plot = "eigplot_" + solution + "_n_" + str(n) +  ".PDF"


os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output.exe " + cla)

plot_path = "/".join([".", "plots", solution])
data_path = "/".join([".", "results", solution])

#Making folders for results and plots
if not os.path.exists(data_path):
    os.makedirs(data_path)
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


os.system("echo creating plots...")
os.system("python3 plot_results.py " + cla)
os.system("echo done.")


os.system("mv" + " numerical" + filename_data + " " + data_path) #Move data file to results directory.
os.system("mv" + " armadillo" + filename_data + " " + data_path)
os.system("mv" + " " + filename_plot + " " + plot_path) #Move file to correct directory.
os.system("mv" + " " + "times_plot_V0.pdf" + " " + plot_path)
