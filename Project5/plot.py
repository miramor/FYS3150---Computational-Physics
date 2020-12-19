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
        sol_met = line[5].strip()
        prob_type = line[6].strip()

    if prob_type == "VD":
        df = pd.read_csv(filename, index_col=False, names=["S", "I", "R", "dD", "n", "B"], skiprows = 1)
    else:
        df = pd.read_csv(filename, index_col=False, names=["S", "I", "R"], skiprows = 1)
    N = df["S"][0] + df["I"][0] + df["R"][0]
    n = len(df["S"])
    return df, N, t, dt, a, b, c, n, sol_met, prob_type

def equilibrium(a,b,c): #Expected (analytic) equilibrium values
    s = b/a
    i = (1-b/a)/(1+b/c)
    r = b/c * (1-b/a) / (1+b/c)
    return s, i, r

def expectation(b,method, cutoff): #Computes expectation values resulting from Monte Carlo simulation
    B = str(b)
    filename = "./Results/" + f"pop_{B}_{method}.csv"
    df, N, t, dt, a, b, c, n, sol_met, prob_type = read_file(filename)
    n_cutoff = int(cutoff*n)

    S_exp = np.mean(df["S"][n_cutoff:])
    I_exp = np.mean(df["I"][n_cutoff:])
    R_exp = np.mean(df["R"][n_cutoff:])

    S_std = np.std(df["S"][n_cutoff:])
    I_std = np.std(df["I"][n_cutoff:])
    R_std = np.std(df["R"][n_cutoff:])

    s_eq, i_eq, r_eq = equilibrium(a,b,c)
    print(f"a = {a}, b = {b}, c = {c}")
    print(f"s* = {s_eq}\ni* = {i_eq}\nr* = {r_eq}")
    return np.array([S_exp, I_exp, R_exp, S_std, I_std, R_std])/N

def Plot_HealthStatus(b, method): #Plot number of susceptibles, infected and immune resulting from MC or RK4
    plt.clf()
    B = str(b)
    filename = "./Results/" + f"pop_{B}_{method}.csv"
    #Reading parameters
    df, N, t, dt, a, b, c, n, sol_met, prob_type = read_file(filename)
    x = np.linspace(0, t, n)
    plt.clf()

    lw = 1
    alpha = 1
    plt.plot(x, df["S"]/N, 'b', alpha = alpha, lw= lw,  label = "S")
    plt.plot(x, df["I"]/N, 'r', alpha = alpha ,lw= lw, label = "I")
    plt.plot(x, df["R"]/N, 'k', alpha = alpha ,lw= lw, label = "R")


    if prob_type == "VD" and sol_met == "MC":
        plt.plot(x, df["n"]/N, alpha = alpha, lw= lw,  label = "Natural deaths")
        plt.plot(x, df["dD"]/N,  alpha = alpha ,lw= lw, label = "Deaths by disease")
        plt.plot(x, df["B"]/N,  alpha = alpha ,lw= lw, label = "Births")

        #deathsDis  = df["dD"].to_numpy()
        #print(f"Death Disease[%], for b={b}: {100*deathsDis[-1]/N:.4f}")
        #print("****************************************")


    #if method == "MC":
    #    plot_avg(x, df, N)

    plt.title(f"Health status, b={b}", size = titlesize)
    plt.ylabel(f"Population [%]", size = labelsize)
    plt.xlabel("Time", size = labelsize)
    plt.legend()
    plt.grid()
    #plt.savefig("./Plots/pop_" + B + "_" + method +  ".pdf")
    plt.savefig(f"./Plots{prob_type}/pop_{B}_{method}_{prob_type}.pdf")


def find_avg(df):
    n = len(df["S"])
    S = np.zeros(n)
    I = np.zeros(n)
    R = np.zeros(n)

    Stot = 0
    Itot = 0
    Rtot = 0

    for i in range(n):
        Stot += df["S"][i]
        Itot += df["I"][i]
        Rtot += df["R"][i]

        S[i] = Stot/(i+1)
        I[i] = Itot/(i+1)
        R[i] = Rtot/(i+1)

    return S, I, R, n

def Plot_HealthStatus2(b, method): #Plot number of susceptibles, infected and immune resulting from two different runs
    plt.clf()
    B = str(b)
    methods = ["MC", "RK4"]
    filenameMC = "./Results/" + f"pop_{B}_{methods[0]}.csv"
    filenameRK4 = "./Results/" + f"pop_{B}_{methods[1]}.csv"
    #Reading parameters
    dfMC, N, t, dtmc, a, b, c, nmc, sol_metmc, prob_type = read_file(filenameMC)
    dfRK4, N, t, dtrk, a, b, c, nrk4, sol_metrk, prob_type = read_file(filenameRK4)

    xmc = np.linspace(0, t, nmc)
    xrk4 = np.linspace(0, t, nrk4)

    plt.clf()

    lw = 1
    alpha = 1

    plt.plot(xmc, dfMC["S"]/N, "b", alpha = alpha, lw= lw,  label = "Vac: S")
    plt.plot(xmc, dfMC["I"]/N, "c", alpha = alpha ,lw= lw, label = "Vac: I")
    plt.plot(xmc, dfMC["R"]/N, "r", alpha = alpha ,lw= lw, label = "Vac: R")
    plt.plot(xrk4, dfRK4["S"]/N, "m", alpha = alpha, lw= lw,  label = "Standard: S", ls = "dashed")
    plt.plot(xrk4, dfRK4["I"]/N, "k", alpha = alpha ,lw= lw, label = "Standard: I", ls = "dashed")
    plt.plot(xrk4, dfRK4["R"]/N, "g", alpha = alpha ,lw= lw, label = "Standard: R", ls = "dashed")

    if prob_type == "VD" and sol_met == "MC":
        plt.plot(x, df["n"]/N, alpha = alpha, lw= lw,  label = "Natural deaths")
        plt.plot(x, df["dD"]/N,  alpha = alpha ,lw= lw, label = "Deaths by disease")
        plt.plot(x, df["B"]/N,  alpha = alpha ,lw= lw, label = "Births")

    #if method == "MC":
    #    plot_avg(x, df, N)

    plt.title(f"Health status, b={b}", size = titlesize)
    plt.ylabel(f"Population [%]", size = labelsize)
    plt.xlabel("Time", size = labelsize)
    plt.legend()
    plt.grid()
    #plt.savefig("./Plots/pop_" + B + "_" + method +  ".pdf")
    plt.savefig(f"./Plots{prob_type}/pop_{B}_{prob_type}.pdf")


def eps_rel(): #prints out (relative) errror between the numerical values and the (analytic) expected ones
    print("\n(Relative) Error and Standard deviation: ")
    cutoffs = [0.50,0.50,0.80,0.50] #needed for equilibtration time for group A, B, C, D respectively

    for i in range(1,5):
        dfRK4, N, t, dt, a, b, c, n, sol_met, prob_type = read_file("./Results/" + f"pop_{i}_RK4.csv")
        s_eq, i_eq, r_eq = equilibrium(a,b,c)
        S_RK4 = dfRK4["S"][n-1]/N
        I_RK4 = dfRK4["I"][n-1]/N
        R_RK4 = dfRK4["R"][n-1]/N
        print(f"______________________________")
        exp_valuesMC = expectation(i, "MC", cutoffs[i-1])
        #print(f"a = {a}, b = {b}, c = {c}")
        #print(f"s* = {s_eq:.4f}\ni* = {i_eq:.4f}\nr* = {r_eq:.4f}")
        if i == 4: #prevent to divide by zero, thus calculating absolute error
            print("RK4:")
            print(f"S_eq = {S_RK4:.4f}  | eps = {abs((S_RK4-s_eq)):.4E} ")
            print(f"I_eq = {I_RK4:.4f}  | eps = {abs((I_RK4-i_eq)):.4E} ")
            print(f"R_eq = {R_RK4:.4f}  | eps = {abs((R_RK4-r_eq)):.4E} ")
            print("MC: ")
            print(f"<S> = {exp_valuesMC[0]:.4f} | STD(S) = {exp_valuesMC[3]:.4E} | eps = {abs((s_eq-exp_valuesMC[0])):.4E}")
            print(f"<I> = {exp_valuesMC[1]:.4f} | STD(I) = {exp_valuesMC[4]:.4E} | eps = {abs((i_eq-exp_valuesMC[1])):.4E}")
            print(f"<R> = {exp_valuesMC[2]:.4f} | STD(R) = {exp_valuesMC[5]:.4E} | eps = {abs((r_eq-exp_valuesMC[2])):.4E}")

        else: #calculating relative error
            print("RK4:")
            print(f"S_eq = {S_RK4:.4f}  | eps_rel = {abs((S_RK4-s_eq)/s_eq):.4E} ")
            print(f"I_eq = {I_RK4:.4f}  | eps_rel = {abs((I_RK4-i_eq)/i_eq):.4E} ")
            print(f"R_eq = {R_RK4:.4f}  | eps_rel = {abs((R_RK4-r_eq)/r_eq):.4E} ")

            print("MC: ")
            print(f"<S> = {exp_valuesMC[0]:.4f} | STD(S) = {exp_valuesMC[3]:.4E} | eps_rel = {abs((s_eq-exp_valuesMC[0])/s_eq):.4E}")
            print(f"<I> = {exp_valuesMC[1]:.4f} | STD(I) = {exp_valuesMC[4]:.4E} | eps_rel = {abs((i_eq-exp_valuesMC[1])/i_eq):.4E}")
            print(f"<R> = {exp_valuesMC[2]:.4f} | STD(R) = {exp_valuesMC[5]:.4E} | eps_rel = {abs((r_eq-exp_valuesMC[2])/r_eq):.4E}")

        # print(f"eps_rel = {abs((exp_valuesMC[0]-s_eq)/s_eq):.4E} ")
        # print(f"eps_rel = {abs((exp_valuesMC[1]-i_eq)/i_eq):.4E} ")
        # print(f"eps_rel = {abs((exp_valuesMC[2]-r_eq)/r_eq):.4E} ")



if __name__ == "__main__":
    prob_type = "prob_type"

    input_correct = False
    while input_correct == False:
        try:
            plot_input = int(input("Would you like to plot with RK4 and MC seperate or together?\n 1. Seperate \n 2. Together \n"))
            if plot_input == 1 or plot_input == 2:
                input_correct = True
            else:
                print("Please enter 1 or 2")
        except:
            print("Please enter 1 or 2")


    # States of S,I,R etc. (most in use)
    for i in range(1,5):
        if plot_input == 1:
            Plot_HealthStatus(i, "RK4")
            Plot_HealthStatus(i, "MC")

        elif plot_input == 2:
            Plot_HealthStatus2(i,"both")


    input_correct = False
    while input_correct == False:
        try:
            print_data = int(input("Would you like to see the error and expectation values? \n 1. Yes \n 2. No \n"))
            if print_data == 1 or print_data == 2:
                input_correct = True
            else:
                print("Please enter 1 or 2")
        except:
            print("Please enter 1 or 2")

    if print_data == 1:
        eps_rel()

        # print(f"______________________________")
        # exp_valuesMC = expectation(i, "MC")
        # print(f"<S> = {exp_valuesMC[0]:.4f} | STD(S) = {exp_valuesMC[3]:.4f}")
        # print(f"<I> = {exp_valuesMC[1]:.4f} | STD(I) = {exp_valuesMC[4]:.4f}")
        # print(f"<R> = {exp_valuesMC[2]:.4f} | STD(R) = {exp_valuesMC[5]:.4f}")

    #Prints out relative error.
