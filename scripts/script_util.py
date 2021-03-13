
import pandas as pd
from os import path 
import json

def parse_path(file_path):
    no_ext = path.splitext(file_path)[0]
    segments = no_ext.split(path.sep)
    relevant = [seg for seg in segments if "=" in seg]
    kv_pairs = [pair for seg in relevant for pair in seg.split("_")]
    kv_pairs = [pair.split("=") for pair in kv_pairs]

    filtered = {} 
    for i, kv in enumerate(kv_pairs):
        if len(kv) < 2:
            continue

        try:
            v = float(kv[1])
        except ValueError:
            v = kv[1]
        filtered[kv[0]] = v

    return filtered 


def read_jsons(json_files):

    result_ls = []
    for jf in json_files:
        
        with open(jf, "r") as f:
            entries = json.load(f)
        
        path_info = parse_path(jf)

        for entry in entries:
            for k, v in path_info.items():
                entry[k] = v

            result_ls.append(entry)

    return result_ls

 
def tabulate_jsons(json_files):

    result_ls = read_jsons(json_files)
    df = pd.DataFrame(result_ls) 
    return df


NICE_NAMES = {"wald_reject": "Rejection Probability (Wald test)",
              "cmh_reject": "Rejection Probability (CMH test)",
              "pat": "$N$",
              "pA": "$p_A$",
              "pB": "$p_B$",
              "N_A": "$N_A$",
              "N_B": "$N_B$",
              "inc": "$\kappa$",
              "fc": "$\lambda_F$",
              "bc": "$\lambda_B$",
              "pr": "v",
              "design": "Trial Design",
              "blockraropt": "BlockRAROpt",
              "traditional": "One-to-One",
              "rar": "RAR",
              "blockrar": "Blocked RAR",
              "excess_failure_frac": "Excess Failure Fraction",
              "excess_failures": "Excess Failures",
              "analysis": "Analysis",
              "wald": "Wald Test",
              "cmh": "CMH Test",
              "cmh_2s": "Two-sided CMH statistic",
              "first_blocksize": "Size of First Block",
              "blocks": "Number of Blocks"
             }

