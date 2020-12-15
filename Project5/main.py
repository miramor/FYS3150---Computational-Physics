import sys
import os
import matplotlib.pyplot as plt

if not os.path.exists("Results"):
    os.makedirs("Results")

if not os.path.exists("Plotsstd"):
    os.makedirs("Plotsstd")
if not os.path.exists("PlotsVD"):
    os.makedirs("PlotsVD")
if not os.path.exists("PlotsVac"):
    os.makedirs("PlotsVac")

cla = " ".join(sys.argv[1:])
print(cla)
os.system("echo Compiling...")
#os.system("make") #compile codes
os.system("make") #compile codes

os.system("echo Executing...")
os.system("./output.exe " + cla)


#os.system("mv" + " e_hist.csv ./Results/")
os.system("echo Making plots...")
os.system("python3 plot.py " + cla)

os.system("echo Full program finished")
