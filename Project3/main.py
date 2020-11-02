import sys
import os
import matplotlib.pyplot as plt

#Input:  system method t_end h
# example:  systemA EC 1 0.0001

if not os.path.exists("Results"):
    os.makedirs("Results")
if not os.path.exists("Plots"):
    os.makedirs("Plots")

cla = sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4]

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output.exe " + cla)

os.system("echo making plots...")
os.system("python3 plot_results.py " + cla)

os.system("echo Full program finished!!")
