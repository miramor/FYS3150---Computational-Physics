import sys
import os
import matplotlib.pyplot as plt

#Input:  system method t_end h
# example:  systemA EC 1 0.0001
cla = sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4]

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output.exe " + cla)

os.system("echo making plots...")
os.system("python3 plot_results.py " + cla)


"""
eps  = 1e-7, t_end = 10 år

dt      errorEuler         errVV
0.01    -6.84863e-05
0.001
1e-4
1e-5
1e-6
1e-7
"""
