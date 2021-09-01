

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from os import path
import script_util as su
from glob import glob
import argparse

def get_blockraropt_qty(score_dir, qty_str, pA, pB, bc="0.001", pr="1.0", stat="scaled_cmh"):

    qty = []
    relevant_files = glob(path.join(score_dir, f"design=blockraropt/fc=*_bc={bc}_pr={pr}_stat={stat}.tsv"))
    fc = []
    
    for fname in relevant_files:
        fc.append(float(fname.split(path.sep)[-1].split("_")[0].split("=")[-1]))
        tab = pd.read_csv(fname, sep="\t")
        tab.set_index(["pA", "pB"], inplace=True)
        qty.append(tab.loc[(pA, pB), qty_str])
        
    srt_idx = np.argsort(fc)
    fc = np.array(fc)[srt_idx]
    qty = np.array(qty)[srt_idx]
    
    return fc, qty


def get_other_qty(score_dir, design_str, qty_str, pA, pB, bc="0.001"):
    
    qty = []
    fc = []
    glob_query = path.join(score_dir, f"design={design_str}/fc=*_bc={bc}_.tsv")
    relevant_files = glob(glob_query)
    
    for fname in relevant_files:
        fc.append(float(fname.split(path.sep)[-1].split("_")[0].split("=")[-1]))
        tab = pd.read_csv(fname, sep="\t")
        tab.set_index(["pA", "pB"], inplace=True)
        qty.append(tab.loc[(pA, pB), qty_str])
        
    srt_idx = np.argsort(fc)
    fc = np.array(fc)[srt_idx]
    qty = np.array(qty)[srt_idx]
    
    return fc, qty


def plot_frontier(xy_list, h_list, label_list, xlabel, ylabel, xl, xu):

    #print("VAR LIST:")
    #print("\t",var_list)
    colors = ["black", "red", "blue", "orange"] 

    #padding = 0.05
    #xmin = min(min(xy[0]) for xy in xy_list)    
    #xmax = max(max(xy[0]) for xy in xy_list)
    #delta_x = xmax - xmin    
    #ymin = min(min(xy[1]) for xy in xy_list)    
    #ymax = max(max(xy[1]) for xy in xy_list)    
    #delta_y = ymax - ymin    

    # Set xlim and ylim
    #plt.xlim(xmin - padding*delta_x, xmax + padding*delta_x)
    if xl is None:
        xl = 0.6
    if xu is None:
        xu = 0.85
    plt.xlim(xl, xu)
    #plt.ylim(ymin - padding*delta_y, ymax + padding*delta_y)
    plt.ylim(-0.1,0.5)

    # plot the zero line (i.e., equal patient allocation)
    plt.plot([-100, 100],[0,0], color="silver", linestyle="--", linewidth=0.5)

    for i, (xy, xyh, label) in enumerate(zip(xy_list, h_list, label_list)):
        plt.errorbar(xy[0], xy[1], xerr=xyh[0], yerr=xyh[1], fmt="none", ecolor=colors[i], elinewidth=0.1, capsize=1, capthick=0.1)
        plt.plot(xy[0], xy[1], color="black", linestyle="--")
        plt.scatter(xy[0], xy[1], color=colors[i], label=su.NICE_NAMES[label])

    plt.xlabel(su.NICE_NAMES[xlabel])
    plt.ylabel(su.NICE_NAMES[ylabel])
    plt.legend()



if __name__=="__main__":

   
    parser = argparse.ArgumentParser()
    parser.add_argument("pA", help="Success probability for treatment A", type=float) 
    parser.add_argument("pB", help="Success probability for treatment B", type=float) 
    parser.add_argument("score_dir", help="directory containing scores")
    parser.add_argument("output_png", help="output PNG file")
    parser.add_argument("--bc", help="blockRARopt block cost parameter",  nargs="+", type=str, default=["0.001"])
    parser.add_argument("--pr", help="blockRARopt prior strength parameter", type=str, default="1.0")
    parser.add_argument("--stat", help="blockRARopt stat parameter", type=str, default="scaled_cmh")
    parser.add_argument("--power_str", help="string identifier for statistical power quantity", required=False, default="cmh_reject")
    parser.add_argument("--outcomes_str", help="string identifier for patient outcome quantity", required=False, default="nA-nB_norm")
    parser.add_argument("--baselines", help="one or more baseline design identifiers", nargs="*", default=["traditional", "rar", "blockrar"], required=False)
    parser.add_argument("--xl", type=float, default=0.6)
    parser.add_argument("--xu", type=float, default=0.85)

    args = parser.parse_args()

    pa, pb = args.pA, args.pB


    xy_ls = []
    h_ls = []
    label_ls = []


    for bc in args.bc:
        fc, power = get_blockraropt_qty(args.score_dir, args.power_str, pa, pb, bc=bc,
                                        pr=args.pr,
                                        stat=args.stat)
        power_h_str = args.power_str + "_h"
        _, power_h = get_blockraropt_qty(args.score_dir, power_h_str, pa, pb, bc=bc,
                                        pr=args.pr,
                                        stat=args.stat)

        _, outcomes = get_blockraropt_qty(args.score_dir, args.outcomes_str, pa, pb, bc=bc,
                                           pr=args.pr,
                                           stat=args.stat)
        outcomes_h_str = args.outcomes_str + "_h"
        _, outcomes_h = get_blockraropt_qty(args.score_dir, outcomes_h_str, pa, pb, bc=bc,
                                           pr=args.pr,
                                           stat=args.stat)
        xy_ls.append((power, outcomes))
        h_ls.append((power_h, outcomes_h))
        label_ls.append("blockraropt")

    bc = args.bc[0]
    for design in args.baselines:
        
        _, power = get_other_qty(args.score_dir, design, args.power_str, pa, pb, bc=bc)
        _, power_h = get_other_qty(args.score_dir, design, power_h_str, pa, pb, bc=bc)
        _, outcomes = get_other_qty(args.score_dir, design, args.outcomes_str, pa, pb, bc=bc)
        _, outcomes_h = get_other_qty(args.score_dir, design, outcomes_h_str, pa, pb, bc=bc)
        xy_ls.append((power, outcomes))
        h_ls.append((power_h, outcomes_h))
        label_ls.append(design)


    plot_frontier(xy_ls, h_ls, label_ls, args.power_str, args.outcomes_str, args.xl, args.xu) 

    plt.savefig(args.output_png, dpi=200)


