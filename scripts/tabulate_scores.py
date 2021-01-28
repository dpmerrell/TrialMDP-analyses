
import pandas as pd
from os import path
import argparse
import json


def parse_path(file_path):
    no_ext = path.splitext(file_path)[0]
    segments = no_ext.split(path.sep)
    relevant = [seg for seg in segments if "=" in seg]
    kv_pairs = [pair for seg in relevant for pair in seg.split("_")]
    kv_pairs = [pair.split("=") for pair in kv_pairs]
    return {kv[0]:kv[1] for kv in kv_pairs}


def read_jsons(json_files):

    result_ls = []
    for jf in json_files:
        
        with open(jf, "r") as f:
            res = json.load(f)
        
        path_info = parse_path(jf)

        for k, v in path_info.items():
            res[k] = v
            

        result_ls.append(res)

    return result_ls

 
def tabulate_data(json_files):

    result_ls = read_jsons(json_files)
    df = pd.DataFrame(result_ls) 
    return df



if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("json_files", help="list of JSON files to tabulate", nargs="+")
    parser.add_argument("output_file", help="CSV file to output")
    
    args = parser.parse_args()

    print(args.output_file)

    df = tabulate_data(args.json_files)
    df.to_csv(args.output_file, sep="\t")


