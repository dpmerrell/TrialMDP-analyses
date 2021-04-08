
import json
import pandas as pd
import argparse
import script_util as su


def compute_scores(summary_df):

    alpha = 0.05
    summary_df["wald_reject"] = (summary_df["wald_p"].astype(float) < alpha).astype(int)
    summary_df["cmh_reject"] = (summary_df["cmh_p"].astype(float) < alpha).astype(int)

    summary_df["wald_2s"] = summary_df["wald"].astype(float).map(lambda x: x**2)
    summary_df["cmh_2s"] = summary_df["cmh"].astype(float).map(lambda x: x**2)

    # MLE estimates
    N_A = summary_df["final_A0"] + summary_df["final_A1"]
    N_B = summary_df["final_B0"] + summary_df["final_B1"]
    summary_df["pA_mle"] = summary_df["final_A1"] / N_A 
    summary_df["pB_mle"] = summary_df["final_B1"] / N_B 
    summary_df["pA_mle_bias"] = summary_df["pA_mle"] - summary_df["pA"]
    summary_df["pB_mle_bias"] = summary_df["pB_mle"] - summary_df["pB"]

    # Total number of patients
    summary_df["pat"] = N_A + N_B    

    summary_df["nA-nB"] = N_A - N_B

    # Excess failures (recall: by construction, pA >= pB)
    summary_df["A_fraction"] = N_A / summary_df["pat"]
    summary_df["excess_failures"] = (summary_df["pA"] - summary_df["pB"])*N_B
    summary_df["failures"] = summary_df["final_A0"] + summary_df["final_B0"]

    # Interim analyses; early stopping
    summary_df["obf_stopping_point"] = (summary_df["interim_n_patients"] / summary_df["pat"])
    summary_df["obf_stopped_early"] = (summary_df["interim_n_patients"] < summary_df["pat"]).astype(float)
    summary_df["obf_reject"] = (summary_df["interim_n_patients"] != summary_df["pat"] | summary_df["cmh_reject"])

    return summary_df



def aggregate_quantities(summary_df, gp_columns):
  
    agg_gp = summary_df.groupby(gp_columns) 
    agg_df = agg_gp.mean() 

    agg_df["nA-nB_10"] = agg_gp["nA-nB"].quantile(0.1)
    agg_df["nA-nB_90"] = agg_gp["nA-nB"].quantile(0.9)

    return agg_df



if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--summaries", required=True, 
                        help="one or more JSON files containing simulation summaries", 
                        nargs="+")
    parser.add_argument("--output_tsv", required=True, help="path to output TSV file")
    parser.add_argument("--gp_columns", nargs="+", default=["pA", "pB", "pat"]) 

    args = parser.parse_args()

    summary_table = su.tabulate_jsons(args.summaries)

    scored_table = compute_scores(summary_table)
    agg_table = aggregate_quantities(scored_table, args.gp_columns)

    agg_table.to_csv(args.output_tsv, sep="\t", index=True)


