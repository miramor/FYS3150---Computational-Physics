import sys
import os
import matplotlib.pyplot as plt

cla = sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4]

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output.exe " + cla)

os.system("echo making plots...")
os.system("python3 plot_results.py " + cla)
