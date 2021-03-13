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


# TODO: use block_size to correctly compute CMH.
#       block_size = -1 for blocked RAR and blockRARopt designs,
#       and some positive integer for RAR designs.
compute_cmh_p_value <- function(history, block_size){

    # compute Cochran-Mantel-Haenszel test statistic
    numerator_sq <- 0.0
    denom <- 0.0
    for(i in 2:length(history)){
        ct <- history[[i]] - history[[i-1]]
        total <- sum(ct)
        rowsum_1 <- ct[1,1] + ct[1,2]
        colsum_1 <- ct[1,1] + ct[2,1]
        numerator_sq <- numerator_sq + ( ct[1,1] - (rowsum_1*colsum_1)/total )
        denom <- denom + (rowsum_1*(total-rowsum_1)*colsum_1*(total-colsum_1)/(total*total*(total-1)))
    }
 
    if(!(is.na(denom)) & denom != 0.0){
        # (two sided)
        #cmh <- numerator_sq*numerator_sq / denom
        # (one sided)
        cmh <- numerator_sq / sqrt(denom)
    }else{
        # (two sided)
        #cmh <- 0.0
        # (one sided)
        #cmh <- -Inf 
        cmh <- 0.0 
    }
  
    # Compute p-value
    # (two sided) 
    #p_value <- pchisq(cmh, 1, lower.tail=FALSE)
    # (one sided) 
    p_value <- pnorm(cmh, lower.tail=FALSE)

    return(list("cmh" = cmh, "p_value" = p_value))
}

interim_analysis <- function(blocks, cur_idx, block_size){
    
    cmh_results <- compute_cmh_p_value(blocks[1:cur_idx], block_size)

    cur_N <- sum(blocks[[cur_idx]])
    final_N <- sum(blocks[[length(blocks)]])

    # TODO: insert code for computing critical value
    critical_value <- 0.05 

    stopped_early <- (cmh[["p_value"]] <= critical_value)
 
    return(stopped_early)
}


simulate_early_stopping <- function(history, block_size){

    interim_info <- list()
    prev_interim <- 0
    n_analyses <- 0
    for(i in 2:length(blocks)){

        cur_N <- sum(blocks[[i]])
        interim_size <- cur_N - prev_interim

        if(interim_size >= block_size){
            stopped_early <- interim_analysis(blocks, i, block_size)
            n_analyses <- n_analyses + 1
            prev_interim <- cur_N
            if(stopped_early){
                break
            }
        }
    }

    interim_info[["n_analyses"]] <- n_analyses
    interim_info[["n_patients"]] <- prev_interim
    final_table <- blocks[[i]]
    interim_info[["final_A0"]] <- final_table[1,2] 
    interim_info[["final_A1"]] <- final_table[1,1] 
    interim_info[["final_B0"]] <- final_table[2,2] 
    interim_info[["final_B1"]] <- final_table[2,1]

    return interim_info
}


analyze_histories <- function (histories, block_size){

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
        cmh <- compute_cmh_p_value(histories[[i]], block_size)
        summary[["cmh"]] <- cmh[["cmh"]]
        summary[["cmh_p"]] <- cmh[["p_value"]]

        interim_info <- simulate_early_stopping(histories[[i]], block_size)
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
                    make_option("--block_size", type="integer", default=-1),
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

block_size <- opts$block_size

summaries <- analyze_histories(histories, block_size)

write_json(summaries, out_json, auto_unbox=TRUE)

