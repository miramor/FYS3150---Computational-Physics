import matplotlib.pyplot as plt
import numpy as np
import sys
n = sys.argv[1]
V = sys.argv[2]
filename = "_n_" + n + "_" + V + ".txt"
n = int(n)
def plot(filename):
    file = open(filename, 'r')
    eigval = np.asarray(file.readline().split(),dtype=float)
    eigvec = np.zeros((n,n))

    lines = file.readlines()
    i = 0
    for line in lines:
        vals = line.split()
        for j in range(n):
            eigvec[j][i] = float(vals[j])
        i += 1

    sorted = np.argsort(eigval)
    return eigval, eigvec, sorted

num_eigval, num_eigvec, num_sort = plot("numerical" + filename)
arma_eigval, arma_eigvec, arma_sort = plot("armadillo" + filename)

n_plots = 2
rho_max = 15
h = rho_max/(n+1)
x = np.asarray([i*h for i in range(n)])
for i in range(n_plots):
    plt.plot(x,num_eigvec[num_sort[i]], label="Num eigenvalue = %.4f" %(num_eigval[num_sort[i]]))
    plt.plot(x,arma_eigvec[arma_sort[i]], label="Arma eigenvalue = %.4f" %(arma_eigval[arma_sort[i]]))

plt.legend()
plt.title("Eigenvectors for %s" %(V))
plt.savefig("eigenvectors_%s_n_%i_plot.PNG" %(V,n))
plt.show()
