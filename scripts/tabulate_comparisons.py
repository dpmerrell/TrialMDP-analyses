
import pandas as pd
import argparse


def load_dfs(tsv_paths):

    return 


def combine_dfs(dfs):

    return


if __name__=="__main__":

    
    parser = argparse.ArgumentParser()

    parser.add_arg("--score_tsvs", help="TSV files containing scores", nargs="+")
    parser.add_arg("--output_tsv", type=str, help="output TSV file")

    args = parser.parse_args()

    dfs = load_dfs(args.score_tsvs)

    table = combine_dfs(dfs)

    table.to_csv(args.output_tsv, sep="\t", index=True)

