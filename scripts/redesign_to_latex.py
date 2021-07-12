
import script_util as su
import pandas as pd
import argparse


if __name__=="__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("score_tsv")
    parser.add_argument("output_tex")
    parser.add_argument("--columns", required=False, nargs="+", default=["pat","blocks","cmh_reject","nA-nB", "nA-nB_10", "nA-nB_90"])
    args = parser.parse_args()

    df = pd.read_csv(args.score_tsv, sep="\t", index_col=[0,1])
    
    df = df[args.columns]   
    df.columns = [su.to_nice(col) for col in df.columns]

    df.to_latex(args.output_tex, float_format="%.3f", 
                index=True, header=True, escape=False)
