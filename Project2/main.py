import os
import sys
import numpy as np


N = int(float(sys.argv[1]))#int(input("Choose n: "))
solution = sys.argv[2]#"V1" #input("Choose a potential V0, V1 or V2: ")

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")

os.system("./output.exe " + str(N) + " " + solution)

os.system("echo creating plots...")
os.system("python3 plot_results.py " + str(N) + " " + solution)

os.system("echo done.")

#Calculating analytical eigenvalues for V0 and V1
ana_eigvalV0 = np.zeros(N)
ana_eigvalV1 = np.zeros(N)
print("\n")
for i in range(N):
    #ana_eigvalV0[i] = (i+1)*np.pi**2
    #ana_eigvalV0[i] = (i+1)**2*np.pi**2
    ana_eigvalV1[i] = 3.0 + 4*i
#print("Analytic eigenvalues: V0\n", ana_eigvalV0)
#print("Analytic eigenvalues: V1 \n", ana_eigvalV1)
