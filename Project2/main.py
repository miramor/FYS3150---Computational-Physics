import os
import sys


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
