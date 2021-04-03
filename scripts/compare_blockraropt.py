
import script_util as su
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


def extract_identifiers(filenames):
    fds = [parse_path(fname) for fname in filenames]
    identifiers = ["fc={}_bc={}_pr={}_stat={}".format(fd["fc"],fd["bc"],fd["pr"],fd["stat"]) for fd in fds]
    return identifiers 


if __name__=="__main__":

    
    parser = argparse.ArgumentParser()

    parser.add_argument("--score_tsvs", help="TSV files containing scores", nargs="+")
    #parser.add_argument("--identifiers", help="string identifiers of trial designs", nargs="+")
    parser.add_argument("--output_xlsx", type=str, help="output XLSX file")
    parser.add_argument("--score_cols", type=str, nargs="+", default=["pat", "cmh_reject", "obf_reject", "obf_stopping_point", "cmh_2s", "nA-nB", "nA-nB_10", "nA-nB_90", "A_fraction", "failures", "excess_failures", "blocks", "first_blocksize", "utility_cmh" ])

    args = parser.parse_args()

    dfs = load_dfs(args.score_tsvs)
    identifiers = extract_identifiers(args.score_tsvs)

    table = combine_dfs(dfs, identifiers, args.score_cols)

    table.to_excel(args.output_xlsx) #, sep="\t", index=True)

