import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

labelsize = 16
titlesize = 18
legendsize = 14


def read_file(filename):
    #Reading parameters
    with open(filename, 'r') as f:
        line = f.readline().split(',')
        t = float(line[0])
        dt = float(line[1])
        a = float(line[2])
        b = float(line[3])
        c = float(line[4])

    df = pd.read_csv(filename, index_col=False, names=["S", "I", "R"], skiprows = 1)
    N = df["S"][0] + df["I"][0] + df["R"][0]
    dp = len(df["S"])
    return df, N, t, dt, a, b, c, dp

def equilibrium(a,b,c):
    s = b/a
    i = (1-b/a)/(1+b/c)
    r = b/c * (1-b/a) / (1+b/c)
    print(f"a = {a}, b = {b}, c = {c}")
    print(f"s* = {s}\n i* = {i}\n r* = {r}")

def expectation(b,method):
    B = str(b)
    filename = "./Results/" + f"pop_{B}_{method}.csv"
    df, N, t, dt, a, b, c, dp = read_file(filename)
    dp_15 = int(0.15*dp)

    S_exp = np.mean(df["S"][dp_15:])
    I_exp = np.mean(df["I"][dp_15:])
    R_exp = np.mean(df["R"][dp_15:])

    S_std = np.std(df["S"][dp_15:])
    I_std = np.std(df["I"][dp_15:])
    R_std = np.std(df["R"][dp_15:])

    equilibrium(a,b,c)
    return np.array([S_exp, I_exp, R_exp, S_std, I_std, R_std])/N

def Plot_HealthStatus(b, method):
    plt.clf()
    B = str(b)
    filename = "./Results/" + f"pop_{B}_{method}.csv"
    #Reading parameters
    df, N, t, dt, a, b, c, dp = read_file(filename)
    x = np.linspace(0, t, dp)
    plt.clf()

    lw = 1
    alpha = 1
    #print(f"{method}: {lw}")
    plt.plot(x, df["S"]/N, 'b', alpha = alpha, lw= lw,  label = "S")
    plt.plot(x, df["I"]/N, 'r', alpha = alpha ,lw= lw, label = "I")
    plt.plot(x, df["R"]/N, 'k', alpha = alpha ,lw= lw, label = "R")


    #if method == "MC":
    #    plot_avg(x, df, N)

    plt.title(f"Health status, b={b}", size = titlesize)
    plt.ylabel(f"Population [%]", size = labelsize)
    plt.xlabel("Time", size = labelsize)
    plt.legend()
    plt.grid()
    plt.savefig("./Plots/pop_" + B + "_" + method + ".pdf")


def plot_avg(x, df, N):
    S, I, R, dp = find_avg(df)
    lw = 1.5
    plt.plot(x, S/N, "b", lw=lw, label= "<S>", ls = "dashed")
    plt.plot(x, I/N, "r", lw = lw,label= "<I>", ls = "dashed")
    plt.plot(x, R/N, "k" , lw = lw, label= "<R>", ls = "dashed")

def find_avg(df):
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

    return S, I, R, dp

def plot_hist(b_val, method):
    plt.clf()
    filename = "./Results/" + f"pop_{b_val}_{method}.csv"
    df, N, t, dt, a, b, c, dp = read_file(filename)
    popGrouped = df.groupby(df["S"],as_index=False).size()
    dp_15 = int(dp*0.15)
    fig, axs = plt.subplots(3)
    axs[0].plot(df["S"][dp_15:])
    axs[1].plot(df["I"][dp_15:])
    axs[3].plot(df["R"][dp_15:])

    sb.distplot(df["S"][dp_15:], norm_hist=True, kde = False, bins = 109)
    #print(f"NORM: {norm.fit(df["E"])})
    plt.title(f"Probability distribution, B={b_val}", size = titlesize)
    plt.ylabel("P(E)", size = labelsize)
    plt.xlabel("Population", size = labelsize)
    plt.savefig("./Plots/hist.pdf")
    print(popGrouped)

#plot_hist(1,"MC")
#for i in [1, 2, 3, 4]:
for i in range(1,5):
    #Plot_HealthStatus(i, "RK4")
    Plot_HealthStatus(i, "MC")
    print(f"______________")
    exp_values = expectation(i, "MC")
    print(f"<S> = {exp_values[0]:.4f} | STD(S) = {exp_values[3]:.4f}")
    print(f"<I> = {exp_values[1]:.4f} | STD(I) = {exp_values[4]:.4f}")
    print(f"<R> = {exp_values[2]:.4f} | STD(R) = {exp_values[5]:.4f}")

#Plot_HealthStatus(2, "MC")




# RK4(2)
# RK4(3)
# RK4(4)
