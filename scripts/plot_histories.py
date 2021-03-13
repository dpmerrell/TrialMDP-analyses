
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import script_util as su
import numpy as np
import matplotlib
import json
import sys


def get_n_patients(history):
    n_patients = [sum([row[0]+row[1] for row in stage]) for stage in history]
    return n_patients

def a_alloc_safe_div(num,denom):
    if denom != 0:
        return num/denom
    elif num == 0 and denom == 0:
        return 0.5
    else:
        raise ValueError

def get_a_alloc(history):
    n_patients = get_n_patients(history)
    a_patients = [stage[0][0]+stage[0][1] for stage in history]

    return [a_alloc_safe_div(a,total) for (a, total) in zip(a_patients, n_patients)]


def get_progress(history):
    n_patients = get_n_patients(history)
    N = n_patients[-1]
    return [pat / N for pat in n_patients]


def plot_history(history):

    progress = get_progress(history) 
    a_alloc = get_a_alloc(history)

    plt.plot(progress, a_alloc, linewidth="0.02", color="gray")

    #plt.scatter(progress, a_alloc)

    return a_alloc[-1] 


def plot_histories(histories, design):

    fig = plt.figure(figsize=(5,3))
    gs = gridspec.GridSpec(nrows=1, ncols=2, 
                           #height_ratios=[1,20], 
                           width_ratios=[10,1])
    fig.subplots_adjust(wspace=0.1, hspace=0.0)

    N = get_n_patients(histories[0])[-1]

    suptitle = "Simulated Histories"
    suptitle += "\n{}; {}={}".format(su.NICE_NAMES[design],
                                     su.NICE_NAMES["pat"], N)
    plt.suptitle(suptitle)
    
    ax1 = fig.add_subplot(gs[0,0])
    plt.xlabel("Trial Progress (fraction of patients)")
    plt.ylabel("Fraction Allocated to Treatment $A$")
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xticks([0,1],["0","1"])
    plt.yticks([0,0.5,1],["0","$\\frac{1}{2}$", "1"])

    final_allocs = []
    for i, h in enumerate(histories):
        #if i % 100 == 0:
        #    print(i, "histories")
        #    print(h)
        #    input()

        final_alloc = plot_history(h)
        final_allocs.append(final_alloc)

    ax2 = fig.add_subplot(gs[0,1])
    plt.hist(final_allocs, bins=25, range=[0,1], orientation="horizontal", color="gray")
    mu = np.mean(final_allocs)
    plt.ylim(0,1)
    plt.xticks([],[])
    ax2.yaxis.tick_right()
    plt.yticks([mu],["{:.2f}".format(mu)])

    return


if __name__=="__main__":

    in_json = sys.argv[1]
    design = sys.argv[2]
    out_png = sys.argv[3]

    with open(in_json, "r") as f:
        d = json.load(f)
        histories = d["histories"]

    plot_histories(histories, design)

    plt.tight_layout()

    plt.savefig(out_png, dpi=200)
        
