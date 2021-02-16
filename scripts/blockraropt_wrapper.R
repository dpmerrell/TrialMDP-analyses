library("blockRARopt")
library("optparse")

# Build the parser
option_list <- list(
    make_option("--alpha_a", default=0.5),
    make_option("--beta_a", default=0.5),
    make_option("--alpha_b", default=0.5),
    make_option("--beta_b", default=0.5),
    make_option("--transition_dist", type="character", default="beta_binom"),
    make_option("--test_statistic", type="character", default="wald")
)

parser <- OptionParser(usage="blockraropt_wrapper.R N_PATIENTS BLOCK_INCREMENT FAILURE_COST BLOCK_COST SQLITE_DB [options]", option_list=option_list)

# Parse the arguments
arguments <- parse_args(parser, positional_arguments=5)
opt <- arguments$options
pargs <- arguments$args


# Run the blockRARopt algorithm
blockRARopt::block_rar_opt(as.integer(pargs[1]),
                           as.integer(pargs[2]), 
                           as.numeric(pargs[3]),
                           as.numeric(pargs[4]), 
                           pargs[5],
                           opt$beta_a, opt$alpha_a,
                           opt$beta_b, opt$alpha_b,
                           transition_dist=opt$transition_dist,
                           test_statistic=opt$test_statistic)

 
