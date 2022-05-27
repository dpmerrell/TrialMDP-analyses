
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import script_util as su
import numpy as np
import matplotlib
import json
import sys


def get_n_patients(history):
    n_patients = [sum(sum(row) for row in stage) for stage in history]
    return n_patients


def get_progress(history):
    n_patients = get_n_patients(history)
    N = n_patients[-1]
    return [pat / N for pat in n_patients]


def get_k(history):
    return len(history) - 1


def plot_history(history):

    k = get_k(history)
    progress = get_progress(history) 

    plt.plot(range(k+1), progress, linewidth=0.02, color="gray")

    return k 


def plot_histories(histories, design):

    fig = plt.figure(figsize=(5,3))
    gs = gridspec.GridSpec(nrows=2, ncols=1, 
                           height_ratios=[2,10]) 
    #fig.subplots_adjust(wspace=0.0, hspace=0.0)

    N = get_n_patients(histories[0])[-1]

    suptitle = "Simulated Histories"
    suptitle += "\n{}; {}={}".format(su.NICE_NAMES[design],
                                     su.NICE_NAMES["pat"], N)
    plt.suptitle(suptitle)
    
    ax1 = fig.add_subplot(gs[1,0])
    plt.xlabel("Trial Block")
    plt.ylabel("Fraction Treated")
    plt.ylim([-0.1,1.1])
    plt.yticks([0,0.5,1],["0","$\\frac{1}{2}$", "1"])

    ks = []
    Kmax = 0
    for i, h in enumerate(histories):
        k = plot_history(h)
        ks.append(k)
        Kmax = max(k, Kmax)
   
    xticks = [0]
    if design in ("blockrar", "blockraropt"):
        for k in range(1,Kmax):
            plt.plot([k,k],[-0.1,1.1], linewidth=0.2, linestyle="--", color="k")
            xticks.append(k)
    xticks.append(Kmax)

    plt.xlim([0,Kmax])
    plt.xticks(xticks, [str(t) for t in xticks])
    plt.plot([0, Kmax], [1.0, 1.0], linewidth=0.5, linestyle="--", color="k")
    plt.plot([0, Kmax], [0.0, 0.0], linewidth=0.5, linestyle="--", color="k")

    ax2 = fig.add_subplot(gs[0,0])
    plt.hist(ks, bins=25, range=[0,Kmax], orientation="vertical", color="gray")
    mu = np.mean(ks)
    plt.xlim([0, Kmax])
    plt.yticks([],[])
    ax2.xaxis.tick_top()
    plt.xticks([mu],["{:.2f}".format(mu)])

    plt.tight_layout()
    
    fig.subplots_adjust(top=0.8)

    return


if __name__=="__main__":

    in_json = sys.argv[1]
    design = sys.argv[2]
    out_png = sys.argv[3]

    with open(in_json, "r") as f:
        d = json.load(f)
        histories = d["histories"]

    plot_histories(histories, design)


    plt.savefig(out_png, dpi=300)
        
