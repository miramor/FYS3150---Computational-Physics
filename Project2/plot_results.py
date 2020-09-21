import matplotlib.pyplot as plt
import numpy as np


def plot(filename):
    file = open("results_num4", 'r')

    x = np.asarray(file.readline().split(),dtype=float)
    eigval = np.asarray(file.readline().split(),dtype=float)
    n = len(x)
    eigvec = np.zeros((n,n))

    lines = file.readlines()
    i = 0
    for line in lines:
        vals = line.split()
        for j in range(n):
            eigvec[i][j] = float(vals[j])
        i += 1
