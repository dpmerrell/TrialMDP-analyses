# summarize_histories.R
# (c) 2021-02 David Merrell
# 
# A script for analyzing trial histories
# (the outputs of simulation_wrapper.R)


library("optparse")
library("jsonlite")


compute_wald_p_value <- function(history){

    # compute Wald test statistic
    ct <- history[[length(history)]]
    p_A <- ct[1,1] / sum(ct[1,1:2])
    p_B <- ct[2,1] / sum(ct[2,1:2])
    p <- sum(ct[1:2,1])/sum(ct)
    N_vec <- rowSums(ct)
    if (0 < p & p < 1){
        # (two sided)
        #W <- (p_A - p_B)^2 / (p*(1-p)) * prod(N_vec) / sum(N_vec)
        # (one sided)
        W <- (p_A - p_B) / sqrt(p*(1-p)) * sqrt(prod(N_vec) / sum(N_vec))
    } else{
        # (two sided)
        #W <- 0
        # (one sided)
        #W <- -Inf
        W <- 0.0
    }
   
    # compute p-value 
    # (two sided)
    #p_value <- pchisq(W, 1, lower.tail=FALSE)
    # (one sided)
    p_value <- pnorm(W, lower.tail=FALSE)
    
    return(list("wald" = W, "p_value" = p_value))
}


compute_cmh_p_value <- function(history, stratum_size){

    # compute Cochran-Mantel-Haenszel test statistic
    numerator_sq <- 0.0
    denom <- 0.0
    stratum_start <- 1
    for(i in 2:length(history)){
        ct <- history[[i]] - history[[stratum_start]]
        total <- sum(ct)
        if((total >= stratum_size) | (i==length(history))){
            rowsum_1 <- ct[1,1] + ct[1,2]
            colsum_1 <- ct[1,1] + ct[2,1]
            numerator_sq <- numerator_sq + ( ct[1,1] - (rowsum_1*colsum_1)/total )
            denom <- denom + (rowsum_1*(total-rowsum_1)*colsum_1*(total-colsum_1)/(total*total*(total-1)))
            
            stratum_start <- i
        }
    }
 
    if(!(is.na(denom)) & denom != 0.0){
        # (one sided)
        cmh <- numerator_sq / sqrt(denom)
    }else{
        cmh <- 0.0 
    }
  
    # Compute p-value
    p_value <- pnorm(cmh, lower.tail=FALSE)

    return(list("cmh" = cmh, "p_value" = p_value))
}

# Estimate effect size in a stratified fashion
compute_stratified_effect <- function(history, stratum_size, pa, pb){

    numerator <- 0.0
    denom <- 0.0

    stratum_start <- 1
    for(i in 2:length(history)){
        ct <- history[[i]] - history[[stratum_start]]
        total <- sum(ct)
        if((total >= stratum_size) | (i==length(history))){
            rowsum_1 <- ct[1,1] + ct[1,2]
            rowsum_2 <- ct[2,1] + ct[2,2] 

            delta <- (ct[1,1]/rowsum_1) - (ct[2,1]/rowsum_2)
            w <- rowsum_1*rowsum_2/(rowsum_1 + rowsum_2) 

            numerator <- numerator + (w*delta)
            denom <- denom + w

            stratum_start <- i
        }
    }

    effect <- numerator / denom
 
    return(effect)
}

interim_analysis <- function(blocks, cur_idx, stratum_size, alpha=0.05){
    
    cmh_results <- compute_cmh_p_value(blocks[1:cur_idx], stratum_size)

    cur_N <- sum(blocks[[cur_idx]])
    final_N <- sum(blocks[[length(blocks)]])

    # For interim analysis, we use alpha-spending
    # with alphas given by O'Brien-Fleming 
    t_frac <- cur_N / final_N
    obf_alpha <- 2.0 - 2.0*pnorm(qnorm(1.0 - (0.5*alpha))/sqrt(t_frac))

    stopped_early <- (cmh_results[["p_value"]] <= obf_alpha)
 
    return(stopped_early)
}


simulate_early_stopping <- function(history, stratum_size){

    interim_info <- list()
    prev_interim <- 0
    n_analyses <- 0
    last_idx <- 1
    for(i in 2:length(history)){

        cur_N <- sum(history[[i]])
        interim_size <- cur_N - prev_interim

        if(interim_size >= stratum_size | i == length(history)){
            stopped_early <- interim_analysis(history, i, stratum_size)
            n_analyses <- n_analyses + 1
            prev_interim <- cur_N
            if(stopped_early){
                last_idx <- i
                break
            }
        }
    }

    interim_info[["n_analyses"]] <- n_analyses
    interim_info[["n_patients"]] <- prev_interim
    final_table <- history[[last_idx]]
    interim_info[["final_A0"]] <- final_table[1,2] 
    interim_info[["final_A1"]] <- final_table[1,1] 
    interim_info[["final_B0"]] <- final_table[2,2] 
    interim_info[["final_B1"]] <- final_table[2,1]

    return(interim_info)
}


analyze_histories <- function (histories, stratum_size, pa, pb){

    summaries <- list()

    for( i in 1:length(histories) ){
        summary <- list()

        # Get the number of blocks, N_A, N_B
        n_stages <- length(histories[[i]])
        summary[["blocks"]] <- (n_stages - 1) 

        summary[["first_blocksize"]] <- sum(histories[[i]][[2]])

        # Store the entries of the final contingency table 
        summary[["final_A0"]] <- histories[[i]][[n_stages]][1,2]
        summary[["final_A1"]] <- histories[[i]][[n_stages]][1,1]
        summary[["final_B0"]] <- histories[[i]][[n_stages]][2,2]
        summary[["final_B1"]] <- histories[[i]][[n_stages]][2,1]

        # Compute Wald statistic and p-value
        wald <- compute_wald_p_value(histories[[i]])
        summary[["wald"]] <- wald[["wald"]]
        summary[["wald_p"]] <- wald[["p_value"]]
 
        # Compute Cochran-Mantel-Haenszel statistic and p-value
        cmh <- compute_cmh_p_value(histories[[i]], stratum_size)
        summary[["cmh"]] <- cmh[["cmh"]]
        summary[["cmh_p"]] <- cmh[["p_value"]]

        summary[["effect_estimate"]] <- compute_stratified_effect(histories[[i]], stratum_size, pa, pb)

        interim_info <- simulate_early_stopping(histories[[i]], stratum_size)
        summary[["interim_n_analyses"]] <- interim_info[["n_analyses"]]
        summary[["interim_n_patients"]] <- interim_info[["n_patients"]]
        summary[["interim_final_A0"]] <- interim_info[["final_A0"]]
        summary[["interim_final_A1"]] <- interim_info[["final_A1"]]
        summary[["interim_final_B0"]] <- interim_info[["final_B0"]]
        summary[["interim_final_B1"]] <- interim_info[["final_B1"]]
 
        summaries[[i]] <- summary 
    }

    return(summaries)
}



# Convert the poorly-formatted read_json data
# to a list of lists of matrices (contingency tables).
# For each contingency table, rows are treatments (1=A and 2=B)
# and columns are outcomes (1=success and 2=failure)
to_matrices <- function(histories){

    result <- list()

    for( h_idx in 1:length(histories) ){
        old_history <- histories[[h_idx]]
        new_history <- list()
        
        for( stage_idx in 1:length(old_history) ){
            old_stage <- old_history[[stage_idx]]
            new_history[[stage_idx]] <- matrix(c(old_stage[[1]][[1]],
                                                 old_stage[[2]][[1]],
                                                 old_stage[[1]][[2]],
                                                 old_stage[[2]][[2]]), 
                                               nrow=2)
        }
        
        result[[h_idx]] <- new_history
    }

    return(result)
}

#################################
# PARSE COMMAND LINE ARGUMENTS
#################################

option_list <- list(
                    make_option("--stratum_frac", type="numeric", default=0.0),
                    make_option("--pa", type="numeric"),
                    make_option("--pb", type="numeric")
                   )

parser <- OptionParser(usage="analyze_histories.R HISTORIES_JSON OUTPUT_JSON", 
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=2)

opts <- arguments$options
pargs <- arguments$args
in_json <- pargs[1]
out_json <- pargs[2]

histories <- read_json(in_json)

histories <- histories$histories

histories <- to_matrices(histories)

h1 <- histories[[1]]
N_max <- sum(h1[[length(h1)]])

stratum_frac <- opts$stratum_frac
stratum_size <- stratum_frac * N_max

summaries <- analyze_histories(histories, stratum_size, opts$pa, opts$pb)

write_json(summaries, out_json, auto_unbox=TRUE, digits=8)


