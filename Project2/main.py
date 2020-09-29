import os
import sys
import numpy as np


N = int(sys.argv[1])#int(input("Choose n: "))
solution = sys.argv[2]#"V1" #input("Choose a potential V0, V1 or V2: ")

os.system("echo compiling...")
os.system("make") #compile codes

os.system("echo executing...")
os.system("./output " + str(N) + " " + solution)

os.system("echo creating plots...")
os.system("python3 plot_results.py " + str(N) + " " + solution)

os.system("echo done.")

#Calculating analytical eigenvalues for V0 and V1
ana_eigvalV0 = np.zeros(N)
ana_eigvalV1 = np.zeros(N)
h = 1/(N+1)
d = 2/(h**2)
a = -1/(h**2)
print("\n")
for i in range(N):
    ana_eigvalV0[i] = d + 2*a*np.cos((i+1)*np.pi/(N+1))
    ana_eigvalV1[i] = 3.0 + 4*i
print("Analytic eigenvalues: V0\n", ana_eigvalV0)
print("Analytic eigenvalues: V1 \n", ana_eigvalV1)
