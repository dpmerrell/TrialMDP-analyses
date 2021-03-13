
import matplotlib.pyplot as plt
import script_util as su
import pandas as pd
import numpy as np
import argparse


def get_score(tsv_file, pA, pB, score_name):

    df = pd.read_csv(tsv_file, sep="\t")
    df.set_index(["pA", "pB"], inplace=True)

    return float(df.loc[(pA, pB), score_name])


def get_N(tsv_file, pA, pB):

    df = pd.read_csv(tsv_file, sep="\t")
    df.set_index(["pA", "pB"], inplace=True)
    
    return int(df.loc[(pA, pB), "pat"])


def collect_scores(tsv_files, pA, pB, score_name):

    path_info = [su.parse_path(tsv) for tsv in tsv_files]
    
    fcs = [pi["fc"] for pi in path_info]
    bcs = [pi["bc"] for pi in path_info]
    scores = [get_score(tsv, pA, pB, score_name) for tsv in tsv_files]

    fc_ls = sorted(list(set(fcs)))
    bc_ls = sorted(list(set(bcs)))
    fc_to_idx = {str(fc): i for i, fc in enumerate(fc_ls)}
    bc_to_idx = {str(bc): i for i, bc in enumerate(bc_ls)}

    score_mat = np.empty((len(fc_ls), len(bc_ls)))
    score_mat[:,:] = np.nan

    for fc, bc, score in zip(fcs, bcs, scores):
        score_mat[fc_to_idx[str(fc)], bc_to_idx[str(bc)]] = score

    return fc_ls, bc_ls, score_mat


def plot_scores(fc_ls, bc_ls, score_mat):

    #plt.imshow(score_mat, vmin=0.0, vmax=1.0, origin="lower", cmap="binary")
    plt.imshow(np.transpose(score_mat), origin="lower", cmap="binary")
    plt.xticks(range(len(fc_ls)), fc_ls)
    plt.yticks(range(len(bc_ls)), bc_ls)

    plt.xlabel(su.NICE_NAMES["fc"])
    plt.ylabel(su.NICE_NAMES["bc"])

 
    return


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("design")
    parser.add_argument("pA", type=float)
    parser.add_argument("pB", type=float)
    parser.add_argument("score_name")
    parser.add_argument("out_png")
    parser.add_argument("--score_tsvs", nargs="+")

    args = parser.parse_args()

    fc, bc, scores = collect_scores(args.score_tsvs,
                                    args.pA, args.pB,
                                    args.score_name)

    N = get_N(args.score_tsvs[0], args.pA, args.pB)

    plot_scores(fc, bc, scores)
    
    plt.colorbar()  
    plt.title("{}\n{}; {}={}; {}={}, {}={}".format(su.NICE_NAMES[args.score_name],
                                        su.NICE_NAMES[args.design],
                                        su.NICE_NAMES["pat"], N,
                                        su.NICE_NAMES["pA"], args.pA,
                                        su.NICE_NAMES["pB"], args.pB)
             )
 
    plt.tight_layout()

    plt.savefig(args.out_png)
 
