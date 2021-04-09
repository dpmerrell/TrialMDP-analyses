
import pandas as pd
import argparse


def load_dfs(tsv_paths):

    result = [] 
    for tsv_p in tsv_paths:
        df = pd.read_csv(tsv_p, sep="\t")
        df.set_index(["pA","pB"], inplace=True)
        result.append(df)

    return result 


def combine_dfs(df_ls, identifiers):

    multicolumns = [(col,ident) for (df, ident) in zip(df_ls, identifiers) for col in df.columns]
    multicolumns = sorted(multicolumns)

    combined = pd.DataFrame(index=df_ls[0].index,
                            columns=pd.MultiIndex.from_tuples(multicolumns))

    for (df, ident) in zip(df_ls, identifiers):
        for col in df.columns:
            combined[(col,ident)] = df[col] 

    return combined


def compute_utility(combined_df, identifiers):
    for ident in identifiers:
        combined_df[("cmh_util", ident)] = combined_df[("cmh_2s", ident)]\
                                         - combined_df[("excess_failures", ident)]*combined_df[("fc","blockraropt")]\
                                         - combined_df[("blocks", ident)]*combined_df[("bc","blockraropt")]

    return combined_df


if __name__=="__main__":

    
    parser = argparse.ArgumentParser()

    parser.add_argument("--agg_tsvs", help="TSV files containing scores", nargs="+")
    parser.add_argument("--identifiers", help="string identifiers of trial designs", nargs="+")
    parser.add_argument("--output_xlsx", type=str, help="output excel file")

    args = parser.parse_args()

    dfs = load_dfs(args.agg_tsvs)

    table = combine_dfs(dfs, args.identifiers)

    table = compute_utility(table, args.identifiers)

    #table.to_csv(args.output_tsv, sep="\t", index=True)
    table.to_excel(args.output_xlsx)


