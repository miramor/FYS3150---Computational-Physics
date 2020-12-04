import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

labelsize = 16
titlesize = 18
legendsize = 14


def Plot_HealthStatus(b, method):
    B = str(b)
    method = method.strip()
    resultpath = "./Results/"
    plotpath = "./Plots/"
    #filename = "pop_" + B + "_MC" +".csv"
    filename = f"pop_{B}_{method}.csv"
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

    lw = 1.5
    alpha = 1
    if method == "MC":
        lw = 0.5
        alpha = 0.5
    #print(f"{method}: {lw}")
    plt.plot(x, df["S"]/N, 'b', alpha = alpha, lw= lw,  label = "S")
    plt.plot(x, df["I"]/N, 'r', alpha = alpha ,lw= lw, label = "I")
    plt.plot(x, df["R"]/N, 'k', alpha = alpha ,lw= lw, label = "R")

    if method == "MC":
        plot_avg(x, df, N)

    plt.title(f"Health status, b={b}", size = titlesize)
    plt.ylabel(f"Population [%]", size = labelsize)
    plt.xlabel("Time", size = labelsize)
    plt.legend()
    plt.grid()
    plt.savefig(plotpath + "pop_" + B + "_" + method + ".pdf")


def plot_avg(x, df, N):
    dp = len(df["S"])
    S = np.zeros(dp)
    I = np.zeros(dp)
    R = np.zeros(dp)

    Stot = 0
    Itot = 0
    Rtot = 0

    for i in range(dp):
        Stot += df["S"][i]
        Itot += df["I"][i]
        Rtot += df["R"][i]

        S[i] = Stot/(i+1)
        I[i] = Itot/(i+1)
        R[i] = Rtot/(i+1)
    lw = 2
    plt.plot(x, S/N, "b", lw=lw, label= "<S>", ls = "dashed")
    plt.plot(x, I/N, "r", lw = lw,label= "<I>", ls = "dashed")
    plt.plot(x, R/N, "k" , lw = lw, label= "<R>", ls = "dashed")

#for i in [1, 2, 3, 4]:
for i in [1]:
    Plot_HealthStatus(i, "RK4")
    Plot_HealthStatus(i, "MC")

# RK4(2)
# RK4(3)
# RK4(4)
