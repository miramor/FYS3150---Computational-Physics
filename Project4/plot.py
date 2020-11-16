import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import numpy as np
from matplotlib.widgets import MultiCursor
from scipy.stats import norm




def plot_d():
    with open("./e_hist.csv", 'r') as f:
        line = f.readline().split(',')
        cutoff = float(line[0])
        numCycles = int(line[1])
        T = float(line[2])
        N = int(line[3])

    df = pd.read_csv("./e_hist.csv", index_col=False, names=["E_mean", "M_mean", "numFlips", "E"], skiprows = 1)

    #plt.rcParams['axes.labelsize'] = 16
    #plt.rcParams['axes.titlesize'] = 16s
    #plt.rcParams['legend.fontsize'] = 19

    #axs[0].set_xlabel('distance (m)') Axis labels
    #axs[0].set_ylabel('Damped oscillation')

    fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex = True, figsize=(4,8))
    E_mean = df["E_mean"].to_numpy()
    M_mean = df["M_mean"].to_numpy()
    E = df["E"].to_numpy()
    numFlips = df["numFlips"].to_numpy()
    MC = np.linspace(cutoff*numCycles,numCycles, E.size)

    axs[0].plot(MC, E_mean, label = "<E>")
    #axs[0].plot(MC, E, label = "E", lw=0.7)
    axs[1].plot(MC, M_mean)
    axs[2].plot(MC, numFlips) #nr of flips
    fig.suptitle(f'2D Isingmodel stabilization at T = {T}', fontsize=14)

    #y-labels
    axs[2].set_xlabel("MonteCarlo cycles")
    axs[0].set_ylabel('Energy')
    axs[1].set_ylabel('Magnetization')
    axs[2].set_ylabel('Total Flips')

    axs[0].legend()
    for i in range(3):
        axs[i].grid()

    #plt.tight_layout()
    plt.savefig(f"stabi_4d_{T}.pdf", dpi = 150)

#plot_d()

def plot_hist():
    plt.clf()
    df = pd.read_csv("./e_hist.csv", index_col=False, names=["E_mean", "M_mean", "numFlips", "E"], skiprows = 1)
    energyGrouped = df.groupby(df["E"],as_index=False).size()
    sb.distplot(df["E"], kde = False)
    plt.ylabel("P(E)")
    plt.xlabel("Energy")
    plt.savefig("prob_E_hist.pdf")
    print(energyGrouped)

#plot_hist()

def plot_obs(L):
    L = str(L)
    path = "./ResultsLong_1mill/"
    savepath = "./Plots/"
    filename = "Observables_" + L + ".csv"
    df = pd.read_csv(path + filename, index_col=False, names=["T", "E", "M", "Cv", "chi"], skiprows = 1)
    #print(df)

    plt.clf()
    plt.plot(df["T"], df["E"], 'ro')
    #plt.plot(df["T"], df["E"])
    plt.axvline(x=2.269)
    plt.title("Energy")
    plt.ylabel("E")
    plt.xlabel("Temperature")
    plt.savefig(savepath + "E_" + L + ".pdf")

    plt.clf()
    plt.plot(df["T"], df["Cv"], 'ro')
    #plt.plot(df["T"], df["Cv"])
    plt.axvline(x=2.269)
    plt.title("Specific heat capacity")
    plt.ylabel("Cv")
    plt.xlabel("Temperature")
    plt.savefig(savepath + "Cv_" + L + ".pdf")

    plt.clf()
    plt.plot(df["T"], df["chi"], 'ro')
    plt.axvline(x=2.269)
    plt.title("Susceptibility")
    plt.ylabel("chi")
    plt.xlabel("Temperature")
    plt.savefig(savepath + "chi_" + L + ".pdf")

    plt.clf()
    plt.plot(df["T"], df["M"], 'ro')
    plt.axvline(x=2.269)
    plt.title("Magnetization")
    plt.ylabel("<|M|>")
    plt.xlabel("Temperature")
    plt.savefig(savepath + "M" + L + ".pdf")



plot_obs(10)

"""
for L in [40, 60, 80, 100, 120]:
    plot_obs(L)
"""




#multi = MultiCursor(fig.canvas, (axs[0], axs[1], axs[2]), color='r', lw=1)
#plt.show()

# Four different plots for opg 4d.
#plt.savefig("E_T1_Ran.pdf", dpi = 350)
#plt.savefig("E_T1_Up.pdf", dpi = 350)
#plt.savefig("E_T2_Ran.pdf", dpi = 350)
#plt.savefig("E_T2_Up.pdf", dpi = 350)

#plt.show()
