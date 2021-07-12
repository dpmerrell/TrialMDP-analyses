library("TrialMDP")
library("optparse")

# Build the parser
option_list <- list(
    make_option("--alpha_a", default=0.5),
    make_option("--beta_a", default=0.5),
    make_option("--alpha_b", default=0.5),
    make_option("--beta_b", default=0.5),
    make_option("--transition_dist", type="character", default="beta_binom"),
    make_option("--test_statistic", type="character", default="wald"),
    make_option("--act_l", type="numeric", default=0.2),
    make_option("--act_u", type="numeric", default=0.8),
    make_option("--act_n", type="integer", default=7)
)

parser <- OptionParser(usage="trialmdp_wrapper.R N_PATIENTS FAILURE_COST BLOCK_COST SQLITE_DB MIN_BLOCK_SIZE BLOCK_INCREMENT  [options]", option_list=option_list)

# Parse the arguments
arguments <- parse_args(parser, positional_arguments=6)
opt <- arguments$options
pargs <- arguments$args


# Run the TrialMDP algorithm
TrialMDP::trial_mdp(as.integer(pargs[1]),
                    as.numeric(pargs[2]),
                    as.numeric(pargs[3]), 
                    pargs[4],
                    as.integer(pargs[5]), 
                    as.integer(pargs[6]),
                    opt$beta_a, opt$alpha_a,
                    opt$beta_b, opt$alpha_b,
                    transition_dist=opt$transition_dist,
                    test_statistic=opt$test_statistic,
                    act_l=opt$act_l, act_u=opt$act_u,
                    act_n=opt$act_n)

 
