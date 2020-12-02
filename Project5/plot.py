import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

labelsize = 16
titlesize = 18
legendsize = 14


def RK4(b):
    B = str(b)
    resultpath = "./Results/"
    plotpath = "./Plots/"
    filename = "pop_" + B + ".csv"
    #Reading parameters
    with open(resultpath + filename, 'r') as f:
        line = f.readline().split(',')
        t = float(line[0])
        dt = float(line[1])
        a = float(line[2])
        b = float(line[3])
        c = float(line[4])


    df = pd.read_csv(resultpath+filename, index_col=False, names=["S", "I", "R"], skiprows = 1)
    N = df["S"][0] + df["I"][0] + df["R"][0]

    x = np.linspace(0, t, len(df["S"]))
    plt.clf()
    plt.plot(x, df["S"]/N, 'b', label = "S")
    plt.plot(x, df["I"]/N, 'r', label = "I")
    plt.plot(x, df["R"]/N, 'k', label = "R")
    plt.title(f"Health status, b={b}", size = titlesize)
    plt.ylabel(f"Population [%]", size = labelsize)
    plt.xlabel("Time", size = labelsize)
    plt.legend()
    plt.grid()
    plt.savefig(plotpath + "pop_" + B + ".pdf")

RK4(1)
RK4(2)
RK4(3)
RK4(4)
