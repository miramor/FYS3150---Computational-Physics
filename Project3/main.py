import sys
import os
import matplotlib.pyplot as plt

#Input:  system method t_end h
# example:  systemA EC 1 0.0001

# 1e-10
#
"""
print(sys.argv[0])
print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
"""

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
0.1     0.00015375      1.98601e-05
0.01    0.000109366     4.03945e-10
0.001   5.8872e-05      1.41106e-14
1e-4    1.21466e-05     9.32414e-18
1e-5    1.38649e-06     5.26922e-17
1e-6    1.40248e-07     9.64669e-17
"""
