import os
import sys
import numpy as np


#n = int(float(sys.argv[1]))#int(input("Choose n: "))
n = int(input("Choose n value:  "))
solution = input("Choose potential type (V0, V1 or V2):  ")
omega = "0"

if solution == "V0":
    rho_max = "1"

if solution == "V1":
    rho_max = "4.6"

if solution == "V2":
    rho_max = "9"
    dict = {"1" : "0.25", "2" : str(1/20), "3" : str(1/54.7386)}
    omegaChoice = input("Choose omega value( 1 = 1/4, 2 = 1/20, 3 = 1/54.7386 ):  ")
    omega = dict[omegaChoice]
    print(omega)

#solution = sys.argv[2]#"V1" #input("Choose a potential V0, V1 or V2: ")
#rho_max = sys.argv[3]
#omega = sys.argv[4]
cla = str(n) + " " + solution + " " + rho_max + " " + omega

if solution == "V2":
    filename_data = "_n_" + str(n) + "_" + solution + "_w_%.3f.txt" %(float(omega))
    #filename_plot1 = "eigplot_" + solution + "_n_" + str(n) + "_w_%.6f.PNG" %(float(omega))
    filename_plot2 = "eigplot_" + solution + "_n_" + str(n) + "_w_%.3f.PDF" %(float(omega))
else:
    filename_data = "_n_" + str(n) + "_" + solution + ".txt"
    #filename_plot1 = "eigplots_" + solution + "_n_" + str(n) +  ".PNG"
    filename_plot2 = "eigplots_" + solution + "_n_" + str(n) +  ".PDF"


os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output.exe " + cla)

plot_path = "/".join([".", "plots", solution])
data_path = "/".join([".", "results", solution])

if not os.path.exists(data_path):
    os.makedirs(data_path)
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


os.system("echo creating plots...")
os.system("python3 plot_results.py " + cla)
os.system("echo done.")

try:
    os.system("mv" + " numerical" + filename_data + " " + data_path) #Move data file to results directory.
    os.system("mv" + " armadillo" + filename_data + " " + data_path)
    os.system("mv" + " " + "TimeTable.csv" + " " + data_path)
    os.system("mv" + " " + filename_plot + " " + plot_path) #Move file to correct directory.
    os.system("mv" + " " + "times_plot_V0.pdf" + " " + plot_path)
except:
    pass
