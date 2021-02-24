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
        W <- (p_A - p_B)^2 / (p*(1-p)) * prod(N_vec) / sum(N_vec)
    } else{
        W <- 0
    }
   
    # compute p-value 
    p_value <- pchisq(W, 1, lower.tail=FALSE)
    
    return(list("wald" = W, "p_value" = p_value))
}


compute_cmh_p_value <- function(history){

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
        cmh <- numerator_sq*numerator_sq / denom
    }else{
        cmh <- 0.0
    }
  
    # Compute p-value 
    p_value <- pchisq(cmh, 1, lower.tail=FALSE)

    return(list("cmh" = cmh, "p_value" = p_value))
}


analyze_histories <- function (histories){

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
        cmh <- compute_cmh_p_value(histories[[i]])
        summary[["cmh"]] <- cmh[["cmh"]]
        summary[["cmh_p"]] <- cmh[["p_value"]]
 
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


parser <- OptionParser(usage="analyze_histories.R HISTORIES_JSON OUTPUT_JSON")#, option_list=option_list)
arguments <- parse_args(parser, positional_arguments=2)

pargs <- arguments$args
in_json <- pargs[1]
out_json <- pargs[2]

histories <- read_json(in_json)
histories <- histories$histories
histories <- to_matrices(histories)

summaries <- analyze_histories(histories)

write_json(summaries, out_json, auto_unbox=TRUE)


