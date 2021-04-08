
import pandas as pd
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
            combined[(col,ident)] = df[col] 

    return combined


if __name__=="__main__":

    
    parser = argparse.ArgumentParser()

    parser.add_argument("--score_tsvs", help="TSV files containing scores", nargs="+")
    parser.add_argument("--identifiers", help="string identifiers of trial designs", nargs="+")
    parser.add_argument("--output_tsv", type=str, help="output TSV file")
    #parser.add_argument("--score_cols", type=str, nargs="+", default=["pat","wald_reject","cmh_reject","wald_2s", "cmh_2s", "pA_mle_bias", "pB_mle_bias", "nA-nB", "nA-nB_10", "nA-nB_90", "A_fraction", "excess_failures", "blocks", "first_blocksize", "utility_wald","utility_cmh" ])
    parser.add_argument("--score_cols", type=str, nargs="+", default=["pat", "cmh_reject", "obf_reject", "obf_stopped_early", "cmh_2s", "nA-nB", "nA-nB_10", "nA-nB_90", "A_fraction", "failures", "excess_failures", "blocks", "first_blocksize", "utility_cmh" ])

    args = parser.parse_args()

    dfs = load_dfs(args.score_tsvs)

    table = combine_dfs(dfs, args.identifiers, args.score_cols)

    table.to_excel(args.output_tsv) #, sep="\t", index=True)

