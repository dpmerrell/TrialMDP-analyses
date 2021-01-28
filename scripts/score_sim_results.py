
import argparse
import json


def type_1_error(sim_results, true_p_a, true_p_b):

    if None in (true_p_a, true_p_b):
        return None
    elif true_p_a != true_p_b:
        return None
    elif true_p_a == true_p_b:
        return sim_results["n_rejected"] / sim_results["trial_count"]


def statistical_power(sim_results, true_p_a, true_p_b):

    if None in (true_p_a, true_p_b):
        return None
    elif true_p_a != true_p_b:
        return sim_results["n_rejected"] / sim_results["trial_count"]
    elif true_p_a == true_p_b:
        return None

def excess_failure_frac(sim_results, true_p_a, true_p_b):
    
    if None in (true_p_a, true_p_b):
        return None

    if true_p_a > true_p_b:
        p_u = true_p_a
        p_l = true_p_b
        n_l = sim_results["N_B"]
    else:
        p_u = true_p_b
        p_l = true_p_a
        n_l = sim_results["N_A"]

    n_patients = sim_results["N_A"] + sim_results["N_B"]   

    xf = (p_u - p_l) * n_l / n_patients
 
    return xf


def compute_scores(sim_results, args):
    
    true_p_a = args.true_p_a
    true_p_b = args.true_p_b
    failure_cost = args.failure_cost
    block_cost = args.block_cost

    sim_results["type_1_error"] = type_1_error(sim_results,
                                               true_p_a,
                                               true_p_b)

    sim_results["power"] = statistical_power(sim_results,
                                             true_p_a,
                                             true_p_b) 

    sim_results["excess_failure_frac"] = excess_failure_frac(sim_results,
                                                         true_p_a,
                                                         true_p_b)
    return sim_results


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_json")
    parser.add_argument("output_json")
    parser.add_argument("--true_p_a", type=float)
    parser.add_argument("--true_p_b", type=float)
    parser.add_argument("--failure_cost", type=float)
    parser.add_argument("--block_cost", type=float)

    args = parser.parse_args()

    with open(args.input_json, "r") as f_in:
        sim_results = json.load(f_in)

    scores = compute_scores(sim_results, args)
 
    with open(args.output_json, "w") as f_out:
        json.dump(scores, f_out)
