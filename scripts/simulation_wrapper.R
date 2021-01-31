# simulation_wrapper.R
# (c) 2021-01 David Merrell
# 
# A script for running monte carlo simulations
# of trial designs. The results are saved to a 
# JSON file.


library("blockRARopt")
library("optparse")
library("jsonlite")

#source(here::here("simulation_util.R"))
#source(here::here("policies.R"))
#source(here::here("summary_updates.R"))


##############################################
# HELPER FUNCTIONS 
##############################################

# Simulate one block of one trial.
# Given a state and action, return the next state.
simulate_block <- function(cur_state, action, true_p_A, true_p_B){
    n_A <- rbinom(1, action[1], true_p_A)
    n_B <- rbinom(1, action[2], true_p_B)
    cur_state[1,1] <- cur_state[1,1] + n_A 
    cur_state[1,2] <- cur_state[1,2] + (action[1] - n_A)
    cur_state[2,1] <- cur_state[2,1] + n_B 
    cur_state[2,2] <- cur_state[2,2] + (action[2] - n_B)
    return(cur_state)
}


# Simulate a whole trial
simulate_trial <- function(N_patients, true_p_A, true_p_B, policy){

    history <- list()
    cur_state <- matrix(rep(0,4), nrow=2)
    cur_block <- 1

    history[[cur_block]] <- cur_state

    while(sum(cur_state) < N_patients){

        action <- policy(cur_state)
        cur_state <- simulate_block(cur_state, action, true_p_A, true_p_B)

        cur_block = cur_block + 1
        history[[cur_block]] <- cur_state

    } 

    return(history)
}


# Simulate N trials
simulate_N_trials <- function(N_trials, N_patients, 
                              true_p_A, true_p_B, 
                              policy, summary_update){

    summary <- list()

    for( trial in 1:N_trials ){
        trial_history <- simulate_trial(N_patients, 
                                        true_p_A, true_p_B, 
                                        policy)
        summary <- summary_update(summary, trial_history) 
    }
    
    return(summary)
}



##############################################
# DEFINE POLICIES
##############################################
traditional_policy <- function(cur_state, N_patients){
    N_A <- round(0.5*N_patients)
    return( c(N_A, N_patients - N_A) )
}



blockraropt_policy <- function(cur_state, conn){
    res = blockRARopt::fetch_result(conn, cur_state[1,2], cur_state[1,1],
                                          cur_state[2,2], cur_state[2,1])
    N <- res[["BlockSize"]]
    N_A <- res[["AAllocation"]]
    return( c(N_A, N - N_A) )
}


##############################################
# DEFINE SUMMARY UPDATES 
##############################################

# After completing a trial, compute some summary information.
# Does a Wald test reject the null hypothesis p_A == p_B?
wald_test_summary_update <- function(summary, history, alpha){

    # Initialize an empty summary
    if(length(summary) == 0){
        summary[["trial_count"]] <- 0
        summary[["n_rejected"]] <- 0
        summary[["N_A"]] <- 0
        summary[["N_B"]] <- 0
        summary[["blocks"]] <- 0
    } 

    # perform Wald test
    ct <- history[[length(history)]]
    p_A <- ct[1,1] / sum(ct[1,1:2])
    p_B <- ct[2,1] / sum(ct[2,1:2])
    p <- sum(ct[1:2,1])/sum(ct)
    N_vec <- rowSums(ct)
    if (0 < p & p < 1){
        W <- (p_A - p_B)^2 / (p*(1-p)) * prod(N_vec) / sum(N_vec)
    } else{
        W <- 0
    }
    # Update summary
   
    tryCatch(
      { 
        if( W >= qchisq(1.0-alpha, 1) ){
            summary[["n_rejected"]] <- summary[["n_rejected"]] + 1
        }
      },
      error = function(condition){
        print("VALUE OF history IS ")
        print(history)
        print("VALUE OF ct IS ", ct)
        print("VALUE OF p IS ", p)
        print("VALUE OF N_vec IS ", N_vec)
        print("VALUE OF W IS ", W)
        stop()
      }
    )
    summary[["N_A"]] <- summary[["N_A"]] + N_vec[1]
    summary[["N_B"]] <- summary[["N_B"]] + N_vec[2]
    summary[["blocks"]] <- summary[["blocks"]] + length(history)-1 
    summary[["trial_count"]] <- summary[["trial_count"]] + 1

    return(summary)
}


##############################################
# RUN SCRIPT 
##############################################

##############################
# Build argument parser
option_list <- list(
    make_option("--design", type="character", default="traditional"),
    make_option("--test", type="character", default="wald"),
    make_option("--blockraropt_db", type="character", default=""),
    make_option("--alpha", default=0.05)
)

parser <- OptionParser(usage="simulation_wrapper.R N_PATIENTS N_SIMULATIONS TRUE_P_A TRUE_P_B OUTPUT_JSON [options]", option_list=option_list)


##############################
# Parse the arguments
arguments <- parse_args(parser, positional_arguments=5)
opt <- arguments$options
pargs <- arguments$args

N_patients <- as.integer(pargs[1])
N_simulations <- as.integer(pargs[2])
true_p_A <- as.numeric(pargs[3])
true_p_B <- as.numeric(pargs[4])
output_json <- pargs[5]


##############################
# Build the trial design
if(opt$design == "traditional"){
    policy <- function(ct){ return(traditional_policy(ct, N_patients)) }
} else if (opt$design == "blockraropt"){
    conn <- blockRARopt::connect_to_results(opt$blockraropt_db)
    policy <- function(ct){ return(blockraropt_policy(ct, conn)) }
} else{
    stop( c("design=", opt$design," is not a valid option.") )
}


##############################
# Build the summary_update function
if(opt$test == "wald"){
    summary_update <- function(s, h) { return(wald_test_summary_update(s,h,opt$alpha)) }
} else{
    stop( c("test=", opt$test," is not a valid option.") )
}


##############################
# Perform the simulations!
results <- simulate_N_trials(N_simulations, N_patients, 
                             true_p_A, true_p_B, 
                             policy, summary_update)



##############################
# Output result JSON
jsonlite::write_json(results, output_json, auto_unbox=TRUE)

