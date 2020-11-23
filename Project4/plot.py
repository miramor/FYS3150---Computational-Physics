import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import numpy as np
from matplotlib.widgets import MultiCursor
from scipy.stats import norm
from sklearn.linear_model import LinearRegression

labelsize = 12
titlesize = 14
legendsize = 12

def plot_stabi():
    with open("./e_hist.csv", 'r') as f:
        line = f.readline().split(',')
        cutoff = float(line[0])
        numCycles = int(line[1])
        T = float(line[2])
        N = int(line[3])

    df = pd.read_csv("./e_hist.csv", index_col=False, names=["E_mean", "M_mean", "numFlips", "E"], skiprows = 1)
    fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex = True, figsize=(4,8))
    E_mean = df["E_mean"].to_numpy()
    M_mean = df["M_mean"].to_numpy()
    E = df["E"].to_numpy()
    numFlips = df["numFlips"].to_numpy()
    MC = np.linspace(cutoff*numCycles,numCycles, E.size)

    axs[0].plot(MC, E, label = "E", lw=0.4)
    axs[0].plot(MC, E_mean, label = "<E>")
    axs[1].plot(MC, M_mean)
    axs[2].plot(MC, numFlips) #nr of flips
    fig.suptitle(f'2D Isingmodel stabilization at T = {T}', fontsize=14)

    #y-labels
    axs[2].set_xlabel("MonteCarlo cycles")
    axs[0].set_ylabel('Energy', size = labelsize)
    axs[1].set_ylabel('Magnetization', size = labelsize)
    axs[2].set_ylabel('Total Flips', size = labelsize)

    axs[0].legend(prop={'size': legendsize})
    for i in range(3):
        axs[i].grid()

    #plt.tight_layout()
    plt.savefig(f"stabi_4d_{T}.pdf", dpi = 150)

    #plt.clf()
    #plt.figure(figsize=(10,8))
    #plt.plot(MC[400*4:10000], E[400*4:10000], label= "Energy 5000 samples")
    #plt.savefig("just_e.pdf", dpi = 300)

#plot_d()

def plot_hist():
    plt.clf()
    df = pd.read_csv("./e_hist.csv", index_col=False, names=["E_mean", "M_mean", "numFlips", "E"], skiprows = 1)
    energyGrouped = df.groupby(df["E"],as_index=False).size()
    sb.distplot(df["E"], norm_hist=True, kde = False, bins = 109)
    #print(f"NORM: {norm.fit(df["E"])})
    plt.title("Probability distribution, T: 2.4", size = titlesize)
    plt.ylabel("P(E)", size = labelsize)
    plt.xlabel("Energy", size = labelsize)
    plt.savefig("prob_E_hist.pdf")
    print(energyGrouped)

def plot_hist_nump():
    plt.clf()
    df = pd.read_csv("./e_hist_up_1.csv", index_col=False, names=["E_mean", "M_mean", "numFlips", "E"], skiprows = 1)
    energyGrouped = df.groupby(df["E"],as_index=False).size()
    variance = np.std(df["E"])**2
    print(energyGrouped.size)

    """
    hist2, bins2 = np.histogram(df["E"], bins = 109)
    bins2 = bins2[:-1] + (bins2[1:]-bins2[:-1])/2
    plt.bar(bins2, hist2)
    plt.savefig("barsBasic.pdf",dpi = 300)
    """
    binSize = energyGrouped.size
    binSize = 109
    L = 20
    hist,bins = np.histogram(df["E"]*400, bins = binSize)
    bins = bins[:-1] + (bins[1:]-bins[:-1])/2
    deltax = np.abs(bins[1]-bins[0])
    a = np.sum(hist*deltax)
    hist = hist/a
    plt.bar(bins/(L*L),hist, width = 1/(binSize))
    plt.annotate(f"Var= {variance:.2e}", xy=(-0.8, 0.006), xytext=(-1.1, 0.006), size = 12)
    plt.title("Probability distribution, T=1", size = titlesize)
    plt.ylabel("P(E)", size = labelsize)
    plt.xlabel("Energy", size = labelsize)
    plt.savefig("barplot_up_1.pdf", dpi = 300)

#plot_hist()
#plot_hist_nump()

def plot_obs(L):
    L = str(L)
    path = "./Results/"
    savepath = "./Plots/"
    filename = "Observables_" + L + ".csv"
    df = pd.read_csv(path + filename, index_col=False, names=["T", "E", "M", "Cv", "chi"], skiprows = 1)
    #print(df.sort_values(by=["T"]))
    df = df.sort_values(by=["T"])

    plt.clf()
    plt.plot(df["T"], df["E"], 'ro')
    #plt.axvline(x=2.269)
    plt.title("Energy", size = titlesize)
    plt.ylabel("E", size = labelsize)
    plt.xlabel("Temperature", size = labelsize)
    plt.grid()
    plt.savefig(savepath + "E_" + L + ".pdf")


    approx_coeff=np.polyfit(df["T"].to_numpy(),df["Cv"].to_numpy(),3)
    approx_func=np.poly1d(approx_coeff)

    plt.clf()
    plt.plot(df["T"], df["Cv"], 'ro')
    plt.plot(df["T"], approx_func(df["T"].to_numpy()))
    #plt.axvline(x=2.269)
    plt.title("Specific heat capacity", size = titlesize)
    plt.ylabel("Cv", size = labelsize)
    plt.grid()
    plt.xlabel("Temperature", size = labelsize)
    plt.savefig(savepath + "Cv_" + L + ".pdf")

    plt.clf()
    plt.plot(df["T"], df["chi"], 'ro')
    #plt.axvline(x=2.269)
    plt.title("Susceptibility", size = titlesize)
    plt.ylabel("chi", size = labelsize)
    plt.grid()
    plt.xlabel("Temperature", size = labelsize)
    plt.savefig(savepath + "chi_" + L + ".pdf")

    plt.clf()
    plt.plot(df["T"], df["M"], 'ro')
    #plt.axvline(x=2.269)
    plt.title("Magnetization", size = titlesize)
    plt.ylabel("<|M|>", size = labelsize)
    plt.grid()
    plt.xlabel("Temperature", size = labelsize)
    plt.savefig(savepath + "M" + L + ".pdf")





"""
for L in [40, 60, 80, 100, 120]:
    plot_obs(L)
"""
def plot_all_obs(L):
    plt.clf()
    path = "./Results/"
    savepath = "./Plots_all/"
    dfs = []
    for i in L:
        filename = "Observables_" + str(i) + ".csv"
        df = pd.read_csv(path + filename, index_col=False, names=["T", "E", "M", "Cv", "chi"], skiprows = 1)
        dfs.append(df.sort_values(by=["T"]))



    linestyle = [":", "-.", "--", "-"]
    color = ["r", "g", "b", "y"]
    alpha = 0.8
    lw = 0.8
    #Plot Cvs
    plt.clf()
    for i in range(len(L)):
        plt.plot(dfs[i]["T"], dfs[i]["Cv"],color[i] + linestyle[i],  label="N=" + str(L[i]), alpha=alpha, linewidth = lw)
        plt.plot(dfs[i]["T"], dfs[i]["Cv"],color[i]+"o", markersize=5)

    #plt.axvline(x=2.269)
    plt.title("Specific heat capacity", size = titlesize)
    plt.ylabel("Cv", size = labelsize)
    plt.xlabel("Temperature", size = labelsize)
    plt.grid()
    plt.legend(prop={'size': legendsize})
    plt.savefig(savepath + "Cv_all.pdf")
    plt.clf()

    for i in range(len(L)):
        plt.plot(dfs[i]["T"], dfs[i]["chi"],color[i] + linestyle[i], label="N="+str(L[i]), alpha=alpha, linewidth = lw)
        plt.plot(dfs[i]["T"], dfs[i]["chi"],color[i]+"o", markersize=5)
    #plt.axvline(x=2.269)
    plt.title("Susceptibility", size = titlesize)
    plt.ylabel("chi", size = labelsize)
    plt.xlabel("Temperature", size = labelsize)
    plt.legend(prop={'size': legendsize})
    plt.grid()
    plt.savefig(savepath + "chi_all.pdf")
    plt.clf()


    for i in range(len(L)):
        dfs[i] = dfs[i].apply (pd.to_numeric, errors='coerce')
        plt.plot(dfs[i]["T"], dfs[i]["M"],color[i] + linestyle[i], label="N="+str(L[i]), alpha=alpha, linewidth = lw)
        plt.plot(dfs[i]["T"], dfs[i]["M"],color[i]+"o", markersize=5)
    #plt.axvline(x=2.269)
    plt.title("Magnetization", size = titlesize)
    plt.legend(prop={'size': legendsize})
    plt.ylabel("<|M|>", size = labelsize)
    plt.grid()
    plt.xlabel("Temperature", size = labelsize)
    plt.savefig(savepath + "M_all.pdf")


def T_critical():
    path = "./Results/"
    savepath = "./Plots_all/"
    dfs = []
    L = np.array([40,60,80,100])
    T_c = np.zeros(len(L))
    Tc_fit = np.zeros(len(L))
    for i in range(len(L)):
        filename = "Observables_" + str(L[i]) + ".csv"
        df = pd.read_csv(path + filename, index_col=False, names=["T", "E", "M", "Cv", "chi"], skiprows = 1)
        T_c[i] = df["T"][np.argmax(df["Cv"].to_numpy())]
        approx_coeff = np.polyfit(df["T"].to_numpy(),df["Cv"].to_numpy(),3)
        approx_func=np.poly1d(approx_coeff)
        x = np. linspace(2,2.4, 1000)
        Tc_fit[i] = x[np.argmax(approx_func(x))]




    approx_coeff=np.polyfit(1/L,Tc_fit,1)
    approx_func=np.poly1d(approx_coeff)
    print(T_c)
    print(Tc_fit)
    print(approx_coeff)

    plt.clf()
    plt.plot(1/L,Tc_fit,'o')
    plt.plot(1/L,approx_func(1/L), label="T=%.4fx+%.4f" %(approx_coeff[0], approx_coeff[1]))
    # for i in range(len(x_values)):
    #     plt.plot(L[i],T_c[i],"o")
    plt.legend(prop={'size': legendsize})
    plt.grid()
    plt.title("Critical Temperature", size = titlesize)
    plt.xlabel("1/Lattice Size", size = labelsize)
    plt.ylabel("Temperature", size = labelsize)
    plt.savefig("./Plots/T_c.pdf")


#T_critical()
#plot_all_obs([40,60,80,100])
"""
L = [40,60,80,100]
for i in L:
    plot_obs(i)
"""
#plot_d()
#plot_hist()
plot_all_obs([40,60,80,100])
#multi = MultiCursor(fig.canvas, (axs[0], axs[1], axs[2]), color='r', lw=1)
#plt.show()

# Four different plots for opg 4d.
#plt.savefig("E_T1_Ran.pdf", dpi = 350)
#plt.savefig("E_T1_Up.pdf", dpi = 350)
#plt.savefig("E_T2_Ran.pdf", dpi = 350)
#plt.savefig("E_T2_Up.pdf", dpi = 350)
#23  2.3000 -1.35775  0.451313  2.102940  61.703400
#21  2.2500 -1.46358  0.678521  1.929050
#plt.show()
