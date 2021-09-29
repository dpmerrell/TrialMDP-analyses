# simulation_wrapper.R
# (c) 2021-01 David Merrell
# 
# A script for running monte carlo simulations
# of trial designs. The results are saved to a 
# JSON file.


library("TrialMDP")
library("optparse")
library("jsonlite")



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
    res = TrialMDP::fetch_result(conn, cur_state[1,2], cur_state[1,1],
                                          cur_state[2,2], cur_state[2,1])
    N <- res[["BlockSize"]]
    N_A <- res[["AAllocation"]]
    return( c(N_A, N - N_A) )
}


blockrar_policy <- function(cur_state, blocksize_map){

    N_A <- cur_state[1,1]+cur_state[1,2]
    N_B <- cur_state[2,1]+cur_state[2,2]

    rar_alloc <- 0.5
    # start RAR policy *after* treating 25% of patients 
    if (N_A != 0 & N_B != 0 & (N_A + N_B) >= 0.25*blocksize_map[["N_patients"]]){
        p_A <- cur_state[1,1] / N_A 
        p_B <- cur_state[2,1] / N_B
        if (0 < p_A & 0 < p_B){ 
            sq_p_A <- sqrt(p_A)
            rar_alloc <- sq_p_A / (sq_p_A + sqrt(p_B))
        }
    }

    blocksize <- blocksize_map[[toString(N_A + N_B)]]

    a_A <- rbinom(1, blocksize, rar_alloc)

    return( c(a_A, blocksize - a_A) )
}


# (Helper function for building the blocksize map)
build_blocksize_map <- function(n_patients, n_blocks){
    blocksize_map <- list()
    pat <- 0
    block_idx <- 1
    rate <- n_patients/n_blocks

    while (pat < n_patients){
        blocksize <- round(block_idx*rate) - round((block_idx-1)*rate)
        blocksize_map[[toString(pat)]] <- blocksize
        
        pat <- pat + blocksize
        block_idx <- block_idx + 1
    }

    blocksize_map["N_patients"] <- n_patients

    return(blocksize_map)
}

##############################################
# DEFINE SUMMARY UPDATES 
##############################################


full_history_summary_update <- function(summary,history){

    if(length(summary) == 0){
        summary[["histories"]] <- list()
        summary[["n_histories"]] <- 0 
    }

    summary[["n_histories"]] <- summary[["n_histories"]] + 1
    summary[["histories"]][[summary[["n_histories"]]]] <- history

    return(summary)
}


##############################################
# RUN SCRIPT 
##############################################

##############################
# Build argument parser
option_list <- list(
    make_option("--design", type="character", default="traditional"),
    make_option("--blockraropt_db", type="character", default=""),
    make_option("--n_blocks", type="integer", default=1),
    make_option("--target_t1", default=0.05),
    make_option("--target_power", default=0.8)
)

parser <- OptionParser(usage="simulation_wrapper.R N_SIMULATIONS TRUE_P_A TRUE_P_B OUTPUT_JSON [options]", option_list=option_list)


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
} else if (opt$design == "rar") {
    blocksize_map <- build_blocksize_map(N_patients, N_patients)
    policy <- function(ct){ return(blockrar_policy(ct, blocksize_map)) }
} else if (opt$design == "blockrar") {
    blocksize_map <- build_blocksize_map(N_patients, opt$n_blocks)
    policy <- function(ct){ return(blockrar_policy(ct, blocksize_map)) }
} else if (opt$design == "blockraropt"){
    conn <- TrialMDP::connect_to_results(opt$blockraropt_db)
    policy <- function(ct){ return(blockraropt_policy(ct, conn)) }
} else{
    stop( c("design=", opt$design," is not a valid option.") )
}


##############################
# Perform the simulations!
results <- simulate_N_trials(N_simulations, N_patients, 
                             true_p_A, true_p_B, 
                             policy, full_history_summary_update)



##############################
# Output result JSON
jsonlite::write_json(results, output_json, auto_unbox=TRUE)


