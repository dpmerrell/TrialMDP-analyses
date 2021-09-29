
import pandas as pd
import numpy as np
import argparse


def load_dfs(tsv_paths):

    result = [] 
    for tsv_p in tsv_paths:
        df = pd.read_csv(tsv_p, sep="\t")
        df.set_index(["pA","pB"], inplace=True)
        result.append(df)

    return result 


def combine_dfs(df_ls, identifiers, score_cols):

    multicolumns = [(col, ident) for col in score_cols for ident in identifiers]

    combined = pd.DataFrame(index=df_ls[0].index,
                            columns=pd.MultiIndex.from_tuples(multicolumns))

    for col in score_cols:
        for (ident, df) in zip(identifiers, df_ls):
            print("\t", ident)
            combined[(col,ident)] = df[col] 

    return combined


def compute_z_scores(df, identifiers, z_cols, baseline="traditional"):

    for zcol in z_cols:
        for method in identifiers:
            df[(zcol+"_z", method)] = (df[(zcol,method)] - df[(zcol,baseline)])/np.sqrt( df[(zcol+"_var",method)]/df[("pat",method)] + df[(zcol+"_var",baseline)]/df[("pat",baseline)] )

    return df


if __name__=="__main__":

    
    parser = argparse.ArgumentParser()

    parser.add_argument("--score_tsvs", help="TSV files containing scores", nargs="+")
    parser.add_argument("--identifiers", help="string identifiers of trial designs", nargs="+")
    parser.add_argument("--output_tsv", type=str, help="output TSV file")
    parser.add_argument("--score_cols", type=str, nargs="+", default=["pat", "cmh_reject", "cmh_reject_h", "obf_reject", "obf_stopped_early", "cmh_2s", "nA-nB_norm", "nA-nB_norm_h", "nA-nB", "nA-nB_h", "nA-nB_05", "nA-nB_95", "effect_estimate", "effect_bias", "effect_bias_h", "A_fraction", "failures", "excess_failures", "blocks", "first_blocksize", "utility_cmh"])
    parser.add_argument("--z_cols", type=str, nargs="+", default=["cmh_reject", "excess_failures", "nA-nB_norm", "utility_cmh"])
    args = parser.parse_args()

    score_cols = args.score_cols + [col+"_var" for col in args.z_cols]
    
    dfs = load_dfs(args.score_tsvs)

    table = combine_dfs(dfs, args.identifiers, score_cols)

    table = compute_z_scores(table, args.identifiers, args.z_cols)

    table.to_excel(args.output_tsv) #, sep="\t", index=True)

