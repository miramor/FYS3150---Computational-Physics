import sys
import os
import matplotlib.pyplot as plt

"""
if not os.path.exists("Results"):
    os.makedirs("Results")
"""
if not os.path.exists("Plots"):
    os.makedirs("Plots")

cla = " ".join(sys.argv[1:])
print(cla)
os.system("echo Compiling...")
os.system("make all") #compile codes

os.system("echo Executing...")
os.system("./output.exe " + cla)

#os.system("echo Making plots...")
#os.system("python3 plot.py " + cla)

os.system("echo Full program finished")
