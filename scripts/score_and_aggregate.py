
import json
import pandas as pd
import argparse
import script_util as su


def compute_scores(summary_df, fc, bc):

    # Reject the null for statistical tests
    # (NOTE: in the future we will also permit 
    #        interim analyses. This will require
    #        access to the full history, and will
    #        be located in the 
    #        "summarize_histories.R" script.)
    alpha = 0.05
    summary_df["wald_reject"] = (summary_df["wald_p"] < alpha).astype(int)
    summary_df["cmh_reject"] = (summary_df["cmh_p"] < alpha).astype(int)

    # MLE estimates
    N_A = summary_df["final_A0"] + summary_df["final_A1"]
    N_B = summary_df["final_B0"] + summary_df["final_B1"]
    summary_df["pA_mle"] = summary_df["final_A1"] / N_A 
    summary_df["pB_mle"] = summary_df["final_B1"] / N_B 

    # Total number of patients
    summary_df["pat"] = N_A + N_B    

    # Excess failures (recall: by construction, pA >= pB)
    summary_df["excess_failures"] = (summary_df["pA"] - summary_df["pB"])*N_B

    # Total utility
    summary_df["utility_wald"] = summary_df["wald"] \
                                 - summary_df["excess_failures"]*fc \
                                 - summary_df["blocks"]*bc
    summary_df["utility_cmh"] = summary_df["cmh"] \
                                - summary_df["excess_failures"]*fc \
                                - summary_df["blocks"]*bc

    return summary_df



def aggregate_quantities(summary_df, gp_columns):
   
    agg_df = summary_df.groupby(gp_columns).mean() 

    return agg_df



if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--summaries", required=True, 
                        help="one or more JSON files containing simulation summaries", 
                        nargs="+")
    parser.add_argument("--output_tsv", required=True, help="path to output TSV file")
    parser.add_argument("--fc", required=True, type=float, help="Failure cost")
    parser.add_argument("--bc", required=True, type=float, help="Block cost")
    parser.add_argument("--gp_columns", nargs="+", default=["pA", "pB", "pat"]) 

    args = parser.parse_args()

    summary_table = su.tabulate_jsons(args.summaries)

    scored_table = compute_scores(summary_table, args.fc, args.bc)
    agg_table = aggregate_quantities(scored_table, args.gp_columns)

    agg_table.to_csv(args.output_tsv, sep="\t", index=True)


