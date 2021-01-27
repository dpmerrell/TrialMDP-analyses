

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


# After completing a trial, compute some summary information.
# Does a Wald test reject the null hypothesis p_A == p_B?
wald_test_summary_update <- function(summary, history){

    # set alpha for the test
    alpha <- 0.05

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
    W <- (p_A - p_B)^2 / (p*(1-p)) * prod(N_vec) / sum(N_vec)
    
    # Update summary
    
    if( W >= qchisq(1.0-alpha, 1) ){
        summary[["n_rejected"]] <- summary[["n_rejected"]] + 1
    }

    summary[["N_A"]] <- summary[["N_A"]] + N_vec[1]
    summary[["N_B"]] <- summary[["N_B"]] + N_vec[2]
    summary[["blocks"]] <- summary[["blocks"]] + length(history)-1 
    summary[["trial_count"]] <- summary[["trial_count"]] + 1

    return(summary)
}


##########################
# Run a simple analysis
##########################

N_patients <- 64
N_simulations <- 10000

true_p_A <- 0.75
true_p_B <- 0.5

# How does a traditional 1:1, single-block policy fare?
traditional_policy <- function(cur_state){
    return( c(0.5*N_patients, 0.5*N_patients) )
}


traditional_results <- simulate_N_trials(N_simulations, N_patients, 
                             true_p_A, true_p_B, 
                             traditional_policy, 
                             wald_test_summary_update)

print("Results for traditional 1:1, single-block design")
print(traditional_results)


# How about a Blocked RAR policy?
conn <- blockRARopt::connect_to_results("results.sqlite")
block_RAR_policy <- function(cur_state){
    res = blockRARopt::fetch_result(conn, cur_state[1,2], cur_state[1,1],
                                          cur_state[2,2], cur_state[2,1])
    N <- res[["BlockSize"]]
    N_A <- res[["AAllocation"]]
    return( c(N_A, N - N_A) )
}


block_RAR_results <- simulate_N_trials(N_simulations, N_patients, 
                             true_p_A, true_p_B, 
                             block_RAR_policy, 
                             wald_test_summary_update)

print("Results for blocked RAR design")
print(block_RAR_results)

