import sys
import os
import matplotlib.pyplot as plt
"""
print("Choose system:\n
1. Sun, Earth \n
2. Sun, Earth, Jupiter \n
3. All planets + pluto \n
4. Sun Earth Mercury")

system = int(input("Input: "))
"""

cla = sys.argv[1] + " " + sys.argv[2]

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output.exe " + cla)

os.system("echo making plots...")
os.system("python3 plot_results.py " + cla)
