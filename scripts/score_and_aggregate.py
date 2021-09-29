
import json
import pandas as pd
import argparse
import script_util as su
import scipy.stats as st


def compute_scores(summary_df, fc, bc):

    alpha = 0.05
    summary_df["wald_reject"] = (summary_df["wald_p"].astype(float) < alpha).astype(int)
    summary_df["cmh_reject"] = (summary_df["cmh_p"].astype(float) < alpha).astype(int)

    summary_df["wald_2s"] = summary_df["wald"].astype(float).map(lambda x: x**2)
    summary_df["cmh_2s"] = summary_df["cmh"].astype(float).map(lambda x: x**2)

    # MLE estimates
    N_A = summary_df["final_A0"] + summary_df["final_A1"]
    N_B = summary_df["final_B0"] + summary_df["final_B1"]
    #summary_df["pA_mle"] = summary_df["final_A1"] / N_A 
    #summary_df["pB_mle"] = summary_df["final_B1"] / N_B 
    #summary_df["pA_mle_bias"] = summary_df["pA_mle"] - summary_df["pA"]
    #summary_df["pB_mle_bias"] = summary_df["pB_mle"] - summary_df["pB"]
    #summary_df["effect"] = summary_df["pA_mle"] - summary_df["pB_mle"]
    #summary_df["true_effect"] = summary_df["pA"] - summary_df["pB"]
    #summary_df["effect_bias"] = summary_df["effect"] - summary_df["true_effect"] 
    summary_df["effect_estimate"] = summary_df["effect_estimate"].astype(float)
    summary_df["effect_bias"] = summary_df["effect_estimate"].astype(float) - (summary_df["pA"] - summary_df["pB"])

    # Total number of patients
    summary_df["pat"] = N_A + N_B    

    summary_df["nA-nB"] = N_A - N_B

    summary_df["nA-nB_norm"] = summary_df["nA-nB"]/summary_df["pat"]

    # Excess failures (recall: by construction, pA >= pB)
    summary_df["A_fraction"] = N_A / summary_df["pat"]
    summary_df["excess_failures"] = (summary_df["pA"] - summary_df["pB"])*(N_B-N_A)/summary_df["pat"]
    summary_df["failures"] = summary_df["final_A0"] + summary_df["final_B0"]

    # Interim analyses; early stopping
    summary_df["obf_stopping_point"] = (summary_df["interim_n_patients"] / summary_df["pat"])
    summary_df["obf_stopped_early"] = (summary_df["interim_n_patients"] < summary_df["pat"]).astype(float)
    summary_df["obf_reject"] = (summary_df["interim_n_patients"] != summary_df["pat"] | summary_df["cmh_reject"])

    #beta = 0.2
    #chisq_beta = qchisq(1.0-beta)

    # Total utility
    #summary_df["utility_wald"] = summary_df["wald"] \
    summary_df["utility_wald"] = summary_df["wald_2s"]/summary_df["pat"]\
                                 - summary_df["excess_failures"]*fc \
                                 - summary_df["blocks"]*bc
    #summary_df["utility_cmh"] = summary_df["cmh"] \
    summary_df["utility_cmh"] = summary_df["cmh_2s"]/summary_df["pat"]\
                                - summary_df["excess_failures"]*fc \
                                - summary_df["blocks"]*bc

    return summary_df



def aggregate_quantities(summary_df, gp_columns):
  
    var_cols = ["cmh_reject", "excess_failures", "utility_cmh", "nA-nB_norm"]
    q_cols = ["nA-nB"]
    ci_cols = ["nA-nB", "cmh_reject", "effect_bias", "nA-nB_norm"]
    ci_conf = 0.90

    agg_gp = summary_df.groupby(gp_columns) 
    agg_df = agg_gp.mean() 

    for col in var_cols:
        agg_df[col+"_var"] = agg_gp[col].var()

    for col in q_cols:
        agg_df[col+"_10"] = agg_gp[col].quantile(0.1)
        agg_df[col+"_90"] = agg_gp[col].quantile(0.9)
        agg_df[col+"_05"] = agg_gp[col].quantile(0.05)
        agg_df[col+"_95"] = agg_gp[col].quantile(0.95)

    for col in ci_cols:
        agg_df[col+"_n"] = agg_gp[col].count()
        agg_df[col+"_sem"] = agg_gp[col].sem()
        agg_df[col+"_h"] = agg_df[col+"_sem"].values * st.t.ppf(0.5*(1.0 + ci_conf), agg_df[col+"_n"].values)

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


