
import pandas as pd
import argparse
import script_util as su

def select_rows(df, scenarios):

    if scenarios == "null":
        idx = [pair[0] == pair[1] for pair in df.index]
    else:
        idx = [pair[0] != pair[1] for pair in df.index]

    return df.loc[idx,:]



def select_cols(df, col_names, suffix):

    supercols = [pair[0] for pair in df.columns]
   
    relevant_supercols = [] 
    for col in col_names:
        if col+suffix in supercols:
            relevant_supercols.append(col+suffix)
        else:
            relevant_supercols.append(col)

    if suffix == "_z":
        methods = ["rar", "blockrar", "blockraropt"]
    else:
        methods = ["traditional", "rar", "blockrar", "blockraropt"]

    cols = [(col, method) for col in relevant_supercols for method in methods]
    cols = [("pat", "traditional")] + cols

    return df[cols]


def rename_cols(df):
    df.columns = df.columns.set_levels( [[su.to_short(col) for col in level] for level in df.columns.levels] )
    return df 


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_xlsx")
    parser.add_argument("output_tex")
    parser.add_argument("--scenarios", default="null")
    parser.add_argument("--columns", nargs="+", default=["cmh_reject", "cmh_reject_h", "effect_bias", "nA-nB", "nA-nB_h", "blocks", "utility_cmh"])
    parser.add_argument("--suffix", nargs="*", default=[""])

    args = parser.parse_args()

    df = pd.read_excel(args.input_xlsx, engine="openpyxl", header=[0,1], index_col=[0,1])

    df = select_rows(df, args.scenarios)

    if len(args.suffix) == 0:
        suffix = ""
    else:
        suffix = args.suffix[0]

    df = select_cols(df, args.columns, suffix)
    df = rename_cols(df)

    df.to_latex(args.output_tex, float_format="%.3f", 
                                 index=True, header=True, escape=False)
    
