
import script_util as su
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import argparse

def plot_heatmap(x_vals, y_vals, arr, ax1, ax2, score, kv):

    plt.imshow(arr, vmin=0.0, vmax=1.0, origin="lower", cmap="binary")

    title = su.NICE_NAMES[score] + "\n"
    for k,v in kv.items():
        title += su.NICE_NAMES[k] + "="
        if isinstance(v, str):
            title += su.NICE_NAMES[v]
        else:
            title += str(v)
        title += "; "
    title = title[:-2]

    plt.title(title)
    plt.xticks(range(len(x_vals)), x_vals)
    plt.yticks(range(len(y_vals)), y_vals)

    plt.xlabel(su.NICE_NAMES[ax1])
    plt.ylabel(su.NICE_NAMES[ax2])

    for j, x in enumerate(x_vals):
        for i, y in enumerate(y_vals):
            if not np.isnan(arr[i,j]):
                plt.text(j,i, "{:.3f}".format(arr[i,j]), ha="center", va="center", color="orange")

    return


def prepare_data(df, plot_column, kv_pairs, ax1, ax2):

    for k,v in kv_pairs.items():
        df = df.loc[ df[k] == v,:]


    aggregates = df[[ax1, ax2, plot_column]].groupby([ax1,ax2]).mean()

    x_vals = df[ax1].unique()
    x_vals.sort()
    y_vals = df[ax2].unique()
    y_vals.sort()

    arr = np.empty((len(x_vals),len(y_vals)))
    arr[:,:] = np.nan

    for i, x_val in enumerate(x_vals):
        for j, y_val in enumerate(y_vals):
            if (x_val, y_val) in aggregates.index:
                # note: (x,y) in image space == (col, row) in index space
                arr[j,i] = aggregates.loc[(x_val, y_val)]

    return x_vals, y_vals, arr


def assemble_kv_pairs(kv_strs):

    result = {}
    for kvs in kv_strs:
        k_v = kvs.split("=")
        k = k_v[0]
        v = k_v[1]
        if not v.isalpha():
            v = float(v)
        result[k] = v 

    return result


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsv", help="path to TSV file containing data to be plotted")
    parser.add_argument("plot_column", help="name of the column to visualize")
    parser.add_argument("output_png", help="path to output PNG file")
    parser.add_argument("--kv_pairs", help="A list of KEY=VALUE pairs specifying which data to plot", nargs="+")
    parser.add_argument("--ax1", help="column for horizontal axis", default="pA")
    parser.add_argument("--ax2", help="column for vertical axis", default="pB")

    args = parser.parse_args()

    df = pd.read_csv(args.input_tsv, sep="\t", index_col=0)

    kv = assemble_kv_pairs(args.kv_pairs)

    x_vals, y_vals, arr = prepare_data(df, args.plot_column, kv, args.ax1, args.ax2)

    plot_heatmap(x_vals, y_vals, arr, args.ax1, args.ax2, args.plot_column, kv)

    plt.colorbar()
    plt.savefig(args.output_png)

