import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

labelsize = 16
titlesize = 18
legendsize = 14

prob_type = "prob_type"

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

def equilibrium(a,b,c):
    s = b/a
    i = (1-b/a)/(1+b/c)
    r = b/c * (1-b/a) / (1+b/c)
    #print(f"a = {a}, b = {b}, c = {c}")
    # print(f"s* = {s}\n i* = {i}\n r* = {r}")
    return s, i, r

def expectation(b,method, cutoff):

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

def Plot_HealthStatus(b, method):
    plt.clf()
    B = str(b)
    filename = "./Results/" + f"pop_{B}_{method}.csv"
    #Reading parameters
    df, N, t, dt, a, b, c, n, sol_met, prob_type = read_file(filename)
    x = np.linspace(0, t, n)
    plt.clf()

    lw = 1
    alpha = 1
    #print(f"{method}: {lw}")
    plt.plot(x, df["S"]/N, 'b', alpha = alpha, lw= lw,  label = "S")
    plt.plot(x, df["I"]/N, 'r', alpha = alpha ,lw= lw, label = "I")
    plt.plot(x, df["R"]/N, 'k', alpha = alpha ,lw= lw, label = "R")


    if prob_type == "VD" and sol_met == "MC":
        plt.plot(x, df["n"]/N, alpha = alpha, lw= lw,  label = "Natural deaths")
        plt.plot(x, df["dD"]/N,  alpha = alpha ,lw= lw, label = "Deaths by disease")
        plt.plot(x, df["B"]/N,  alpha = alpha ,lw= lw, label = "Births")

        deathsDis  = df["dD"].to_numpy()
        print(f"Death Disease[%], for b={b}: {100*deathsDis[-1]/N:.4f}")
        print("****************************************")


    #if method == "MC":
    #    plot_avg(x, df, N)

    plt.title(f"Health status, b={b}", size = titlesize)
    plt.ylabel(f"Population [%]", size = labelsize)
    plt.xlabel("Time", size = labelsize)
    plt.legend()
    plt.grid()
    #plt.savefig("./Plots/pop_" + B + "_" + method +  ".pdf")
    plt.savefig(f"./Plots{prob_type}/pop_{B}_{method}_{prob_type}.pdf")


def plot_avg(x, df, N):
    S, I, R, n = find_avg(df)
    lw = 1.5
    plt.plot(x, S/N, "b", lw=lw, label= "<S>", ls = "dashed")
    plt.plot(x, I/N, "r", lw = lw,label= "<I>", ls = "dashed")
    plt.plot(x, R/N, "k" , lw = lw, label= "<R>", ls = "dashed")

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

def plot_hist(b_val, method):
    plt.clf()
    filename = "./Results/" + f"pop_{b_val}_{method}.csv"
    df, N, t, dt, a, b, c, n, sol_met, prob_type = read_file(filename)
    popGrouped = df.groupby(df["S"],as_index=False).size()
    n_cutoff = int(n*0.15)
    # fig, axs = plt.subplots(3)
    # axs[0].plot(df["S"][n_cutoff:])
    # axs[1].plot(df["I"][n_cutoff:])
    # axs[2].plot(df["R"][n_cutoff:])

    sb.distplot(df["S"][n_cutoff:]/N, norm_hist=True, kde = False, bins = 109)
    #print(f"NORM: {norm.fit(df["E"])})
    plt.title(f"Probability distribution, B={b_val}", size = titlesize)
    plt.ylabel("P(E)", size = labelsize)
    plt.xlabel("Population [%]", size = labelsize)
    plt.savefig("./Plots/hist.pdf")
    print(popGrouped)

def Plot_PhasePortrait(b, method):
    plt.clf()
    B = str(b)
    filename = "./Results/" + f"pop_{B}_{method}.csv"
    df, N, t, dt, a, b, c, n, sol_met, prob_type = read_file(filename)
    plt.clf()
    lw = 1
    alpha = 1
    #print(f"{method}: {lw}")
    plt.plot(df["S"]/N,  df["I"]/N, 'b', alpha = alpha, lw= lw,  label = "S")
    plt.title(f"Phase portrait, b={b}", size = titlesize)
    plt.ylabel("I/N", size = labelsize)
    plt.xlabel("S/N", size = labelsize)
    #plt.legend()
    plt.grid()
    plt.savefig(f"./Plots/phasepor_{B}_{method}_{prob_type}.pdf")

# for i in range(1,5):
#     Plot_PhasePortrait(i, "RK4")
#plot_hist(1,"MC")
#for i in [1, 2, 3, 4]:
for i in range(1,5):
    a = 2 #filler, not used. To avoid error with for func when running
    #Plot_HealthStatus(i, "RK4")
    #Plot_HealthStatus(i, "MC")
    # print(f"______________________________")
    # exp_valuesMC = expectation(i, "MC")
    # print(f"<S> = {exp_valuesMC[0]:.4f} | STD(S) = {exp_valuesMC[3]:.4f}")
    # print(f"<I> = {exp_valuesMC[1]:.4f} | STD(I) = {exp_valuesMC[4]:.4f}")
    # print(f"<R> = {exp_valuesMC[2]:.4f} | STD(R) = {exp_valuesMC[5]:.4f}")


def eps_rel():
    print("\n(Relative) Error and Standard deviation: ")
    cutoffs = [0.50,0.50,0.80,0.50]
    #for i,cut in zip(range(1,5), cutoffs):

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

        else: #calculating relative error
            print("RK4:")
            print(f"S_eq = {S_RK4:.4f}  | eps_rel = {abs((S_RK4-s_eq)/s_eq):.4E} ")
            print(f"I_eq = {I_RK4:.4f}  | eps_rel = {abs((I_RK4-i_eq)/i_eq):.4E} ")
            print(f"R_eq = {R_RK4:.4f}  | eps_rel = {abs((R_RK4-r_eq)/r_eq):.4E} ")

        print("MC: ")
        print(f"<S> = {exp_valuesMC[0]:.4f} | STD(S) = {exp_valuesMC[3]:.4f} | eps_RK4 = {abs((S_RK4-exp_valuesMC[0])/S_RK4):.4E}")
        print(f"<I> = {exp_valuesMC[1]:.4f} | STD(I) = {exp_valuesMC[4]:.4f} | eps_RK4 = {abs((R_RK4-exp_valuesMC[1])/I_RK4):.4E}")
        print(f"<R> = {exp_valuesMC[2]:.4f} | STD(R) = {exp_valuesMC[5]:.4f} | eps_RK4 = {abs((I_RK4-exp_valuesMC[2])/R_RK4):.4E}")
        #HUSK n = 30%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print(f"eps_rel = {abs((exp_valuesMC[0]-s_eq)/s_eq):.4E} ")
        print(f"eps_rel = {abs((exp_valuesMC[1]-i_eq)/i_eq):.4E} ")
        print(f"eps_rel = {abs((exp_valuesMC[2]-r_eq)/r_eq):.4E} ")

def plot_deathInfected(filename, colNames):
    df = pd.read_csv(filename, index_col=False, names= colNames)
    plt.title("Deaths due vs vaccination constant", size = titlesize)
    plt.plot(df["dI"], df["A"]/400, label = "Group A")
    plt.plot(df["dI"], df["B"]/400, label = "Group B")
    plt.plot(df["dI"], df["C"]/400, label = "Group C")
    plt.plot(df["dI"], df["D"]/400, label = "Group D")
    plt.xlabel("f", size = labelsize)
    plt.ylabel("Deaths[%]", size = labelsize)
    plt.ylim([0,1])
    plt.legend()
    plt.grid()
    plt.savefig("deathInfected.pdf", dpi = 300)
    plt.show()

#plot_deathInfected("death_infected.csv")
colnames1 = ["dI", "A", "B", "C", "D"]
colnames2 = ["dI", "A", "B", "C", "D", "S_fin", "I_fin", "R_fin"]
plot_deathInfected("dI_data.csv", colnames1)


#eps_rel()
#plot_hist(1, "MC")
