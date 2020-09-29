import os
import sys
import numpy as np


n = int(float(sys.argv[1]))#int(input("Choose n: "))
solution = sys.argv[2]#"V1" #input("Choose a potential V0, V1 or V2: ")
rho_max = sys.argv[3]
omega = sys.argv[4]
cla = str(n) + " " + solution + " " + rho_max + " " + omega

if solution == "V2":
    filename_data = "_n_" + str(n) + "_" + solution + "_w_%.6f.txt" %(float(omega))
    filename_plot = "eigplot_" + solution + "_n_" + str(n) + "_w_%.6f.PNG" %(float(omega))
else:
    filename_data = "_n_" + str(n) + "_" + solution + ".txt"
    filename_plot = "eigplots_" + solution + "_n_" + str(n) +  ".PNG"


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


#Calculating analytical eigenvalues for V0 and V1
ana_eigvalV0 = np.zeros(n)
ana_eigvalV1 = np.zeros(n)
print("\n")
for i in range(n):
    #ana_eigvalV0[i] = (i+1)*np.pi**2
    #ana_eigvalV0[i] = (i+1)**2*np.pi**2
    ana_eigvalV1[i] = 3.0 + 4*i
#print("Analytic eigenvalues: V0\n", ana_eigvalV0)
#print("Analytic eigenvalues: V1 \n", ana_eigvalV1)
