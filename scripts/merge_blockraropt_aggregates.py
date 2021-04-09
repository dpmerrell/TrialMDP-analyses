
import script_util as su
import pandas as pd
import numpy as np
import argparse
import json

def load_dfs(tsv_paths, params_ls):

    result = [] 
    for (tsv_p, params) in zip(tsv_paths, params_ls):
        df = pd.read_csv(tsv_p, sep="\t")
        df.set_index(["pA","pB"], inplace=True)
        
        for k, v in params.items():
            df[k] = v

        result.append(df)

    return result 


def fnames_to_params(filenames):
    params = [su.parse_path(fname) for fname in filenames]
    return params


def get_null_df_idx(param_ls, null_fc=0.5, null_bc=0.0, null_pr=0.5):
    for i, p in enumerate(param_ls):
        if p["fc"] == null_fc and p["bc"] == null_bc and p["pr"] == null_pr:
            return i
    raise ValueError


def merge_dfs(df_ls, null_df_idx):
    
    best_df = pd.DataFrame(index=df_ls[0].index, columns=df_ls[0].columns)

    for scenario in df_ls[0].index:
        pA = scenario[0]
        pB = scenario[1]
        if pA != pB:
            relevant_idx = [i for i, df in enumerate(df_ls) if df.loc[scenario, "blocks"] > 1]
            relevant_power = [df_ls[idx].loc[scenario, "cmh_reject"] for idx in relevant_idx]
            best_idx = np.argmax(relevant_power)
            best_df.loc[scenario,:] = df_ls[relevant_idx[best_idx]].loc[scenario,:]
        else:
            best_df.loc[scenario,:] = df_ls[null_df_idx].loc[scenario,:]
            
    return best_df


if __name__=="__main__":

    
    parser = argparse.ArgumentParser()
    parser.add_argument("--agg_tsvs", help="TSV files containing aggregates", nargs="+")
    parser.add_argument("--output_tsv", type=str, help="output TSV containing best results for each scenario")
    parser.add_argument("--null_fc", type=float, default=0.5)
    parser.add_argument("--null_bc", type=float, default=0.0)
    parser.add_argument("--null_pr", type=float, default=0.5)
    args = parser.parse_args()

    blockraropt_params = fnames_to_params(args.agg_tsvs)
    null_df_idx = get_null_df_idx(blockraropt_params, args.null_fc,
                                                      args.null_bc,
                                                      args.null_pr)
    blockraropt_dfs = load_dfs(args.agg_tsvs, blockraropt_params)
    
    best_df = merge_dfs(blockraropt_dfs, null_df_idx)
 
    best_df.to_csv(args.output_tsv, sep="\t", index=True)


