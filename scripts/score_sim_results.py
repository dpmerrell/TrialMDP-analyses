
import argparse
import json


## After completing a trial, compute some summary information.
## Does a Wald test reject the null hypothesis p_A == p_B?
#wald_test_summary_update <- function(summary, history, alpha){
#
#    # Initialize an empty summary
#    if(length(summary) == 0){
#        summary[["trial_count"]] <- 0
#        summary[["n_rejected"]] <- 0
#        summary[["N_A"]] <- 0
#        summary[["N_B"]] <- 0
#        summary[["blocks"]] <- 0
#    } 
#
#    # perform Wald test
#    ct <- history[[length(history)]]
#    p_A <- ct[1,1] / sum(ct[1,1:2])
#    p_B <- ct[2,1] / sum(ct[2,1:2])
#    p <- sum(ct[1:2,1])/sum(ct)
#    N_vec <- rowSums(ct)
#    if (0 < p & p < 1){
#        W <- (p_A - p_B)^2 / (p*(1-p)) * prod(N_vec) / sum(N_vec)
#    } else{
#        W <- 0
#    }
#    # Update summary
#   
#    tryCatch(
#      { 
#        if( W >= qchisq(1.0-alpha, 1) ){
#            summary[["n_rejected"]] <- summary[["n_rejected"]] + 1
#        }
#      },
#      error = function(condition){
#        print("VALUE OF history IS ")
#        print(history)
#        print("VALUE OF ct IS ", ct)
#        print("VALUE OF p IS ", p)
#        print("VALUE OF N_vec IS ", N_vec)
#        print("VALUE OF W IS ", W)
#        stop()
#      }
#    )
#    summary[["N_A"]] <- summary[["N_A"]] + N_vec[1]
#    summary[["N_B"]] <- summary[["N_B"]] + N_vec[2]
#    summary[["blocks"]] <- summary[["blocks"]] + length(history)-1 
#    summary[["trial_count"]] <- summary[["trial_count"]] + 1
#
#    return(summary)
#}
#
#
#cmh_test_summary_update <- function(summary, history, alpha){
#
#    # Initialize an empty summary
#    if(length(summary) == 0){
#        summary[["trial_count"]] <- 0
#        summary[["n_rejected"]] <- 0
#        summary[["N_A"]] <- 0
#        summary[["N_B"]] <- 0
#        summary[["blocks"]] <- 0
#    } 
#
#    # perform Cochran-Mantel-Haenszel test
#    numerator_sq = 0.0
#    denom = 0.0
#    for(i in 2:length(history)){
#        ct <- history[[i]] - history[[i-1]]
#        total <- sum(ct)
#        rowsum_1 <- ct[1,1] + ct[1,2]
#        colsum_1 <- ct[1,1] + ct[2,1]
#        numerator_sq <- numerator_sq + ( ct[1,1] - (rowsum_1*colsum_1)/total )
#        denom <- denom + (rowsum_1*(total-rowsum_1)*colsum_1*(total-colsum_1)/(total*total*(total-1)))
#    } 
#    if (denom != 0.0){
#        cmh <- numerator_sq*numerator_sq / denom
#    } else{
#        cmh <- 0.0
#    }
#    # Update summary
#   
#    tryCatch(
#      { 
#        if( cmh >= qchisq(1.0-alpha, 1) ){
#            summary[["n_rejected"]] <- summary[["n_rejected"]] + 1
#        }
#      },
#      error = function(condition){
#        print("VALUE OF history IS ")
#        print(history)
#        print("VALUE OF numerator_sq IS ", numerator_sq)
#        print("VALUE OF denom IS ", denom)
#        print("VALUE OF cmh IS ", cmh)
#        stop()
#      }
#    )
#    N_vec <- rowSums(history[[length(history)]])
#    summary[["N_A"]] <- summary[["N_A"]] + N_vec[1]
#    summary[["N_B"]] <- summary[["N_B"]] + N_vec[2]
#    summary[["blocks"]] <- summary[["blocks"]] + length(history)-1 
#    summary[["trial_count"]] <- summary[["trial_count"]] + 1
#
#    return(summary)
#}



def prob_reject(sim_results):

    return sim_results["n_rejected"] / sim_results["trial_count"]


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

    sim_results["prob_reject"] = prob_reject(sim_results)

    sim_results["excess_failure_frac"] = excess_failure_frac(sim_results,
                                                         true_p_a,
                                                         true_p_b)

    return sim_results


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("histories_json")
    parser.add_argument("output_json")
    parser.add_argument("--true_p_a", type=float)
    parser.add_argument("--true_p_b", type=float)
    parser.add_argument("--failure_cost", type=float)
    parser.add_argument("--block_cost", type=float)

    args = parser.parse_args()

    with open(args.histories_json, "r") as f_in:
        sim_histories = json.load(f_in)

    scores = compute_scores(sim_histories, args)
 
    with open(args.output_json, "w") as f_out:
        json.dump(scores, f_out)
