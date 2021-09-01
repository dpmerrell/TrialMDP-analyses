

from os import path
from scipy.stats import norm

configfile: "config.yaml"

SCRIPT_DIR = "scripts"
RESULT_DIR = config["result_dir"]

####################################################
# SIMULATION STUDY PARAMS
####################################################

#SIM_RESULT_DIR = path.join(RESULT_DIR, "simulation_study")
SIM_RESULT_DIR = RESULT_DIR
POLICY_DB_DIR = path.join(SIM_RESULT_DIR, "precomputed_policies")
HISTORY_DIR = path.join(SIM_RESULT_DIR, "histories")
SUMMARY_DIR = path.join(SIM_RESULT_DIR, "summaries") 
SCORE_DIR = path.join(SIM_RESULT_DIR, "scores")
COMPARISON_DIR = path.join(SIM_RESULT_DIR, "comparisons")
LATEX_COMPARISON_DIR = path.join(SIM_RESULT_DIR, "comparisons_latex")
FIG_DIR = path.join(SIM_RESULT_DIR, "figures")

SIM_PARAMS = config["simulation_params"]
N_SAMPLES = SIM_PARAMS["n_samples"]
TRUE_PA = SIM_PARAMS["true_pa"]
TRUE_PB = SIM_PARAMS["true_pb"]
NULL_HYP_N = SIM_PARAMS["null_hypothesis_n"]
TRUE_P_PAIRS = ["pA={}_pB={}".format(pa, pb) for pa, pb in zip(TRUE_PA,TRUE_PB)]
SIM_ALPHA = SIM_PARAMS["alpha"]
SIM_BETA = SIM_PARAMS["beta"]

BASELINES = config["baseline_designs"]

BLOCKRAROPT_PARAMS = config["trialmdp_params"]
FAILURE_COST = BLOCKRAROPT_PARAMS["failure_cost"]
ACT_N = BLOCKRAROPT_PARAMS["act_n"]
BLOCK_COST = BLOCKRAROPT_PARAMS["block_cost"]
PRIOR_STRENGTH = BLOCKRAROPT_PARAMS["prior_strength"]
OPT_STAT = BLOCKRAROPT_PARAMS["stat"]

BLOCKRAR_PARAMS = config["blockrar_params"]
BLOCKRAR_N_BLOCKS = BLOCKRAR_PARAMS["n_blocks"]


##########################################
# TRIAL REDESIGN PARAMS
##########################################

REDESIGN_DIR = path.join(RESULT_DIR, "trial_redesign")
RED_POLICY_DB_DIR = path.join(REDESIGN_DIR, "precomputed_policies")
RED_HISTORY_DIR = path.join(REDESIGN_DIR, "histories")
RED_SUMMARY_DIR = path.join(REDESIGN_DIR, "summaries") 
RED_SCORE_DIR = path.join(REDESIGN_DIR, "sweep_scores")
RED_TEST_SCORE_DIR = path.join(REDESIGN_DIR, "test_scores")
RED_COMPARISON_DIR = path.join(REDESIGN_DIR, "comparisons")
RED_FIG_DIR = path.join(REDESIGN_DIR, "figures")

REDESIGN_PARAMS = config["trial_redesign_params"]
REDESIGN_N = REDESIGN_PARAMS["n_patients"]
REDESIGN_ALPHA = REDESIGN_PARAMS["alpha"]
REDESIGN_BETA = REDESIGN_PARAMS["beta"]
REDESIGN_PA = REDESIGN_PARAMS["true_pa"]
REDESIGN_PB = REDESIGN_PARAMS["true_pb"]
REDESIGN_SIMS = REDESIGN_PARAMS["n_simulations"]
REDESIGN_P_PAIRS = ["pA={}_pB={}".format(pa, pb) for pa, pb in zip(REDESIGN_PA,REDESIGN_PB)]

REDESIGN_FC = REDESIGN_PARAMS["failure_cost"]
REDESIGN_BC = REDESIGN_PARAMS["block_cost"]
REDESIGN_MINSIZE = REDESIGN_PARAMS["min_blocksize"]
REDESIGN_INCR = REDESIGN_PARAMS["block_incr"]
REDESIGN_STAT = REDESIGN_PARAMS["stat"]
REDESIGN_PR = REDESIGN_PARAMS["prior_strength"]

REDESIGN_CHOSEN_FC = REDESIGN_PARAMS["chosen_fc"]
REDESIGN_CHOSEN_BC = REDESIGN_PARAMS["chosen_bc"]

###############################
# COMPUTE NUMBER OF PATIENTS
def compute_n_patients(pA, pB, alpha=SIM_ALPHA, beta=SIM_BETA):
    
    if pA == pB:
        return NULL_HYP_N

    delta = abs(pA - pB)

    var_b = pB*(1.0-pB)
    var_a = pA*(1.0-pA)

    k = 1.0 

    # one-sided:
    n_a = (norm.ppf(1.0-alpha) + norm.ppf(1.0-beta))**2.0 * (var_a + var_b/k) / (delta*delta)
    # two-sided:
    #n_a = (norm.ppf(1.0 - 0.5*alpha) + norm.ppf(1.0-beta))**2.0 * (var_a + var_b/k) / (delta*delta)
    n_b = k * n_a

    # make it an even number
    n = (int(n_a + n_b) // 2)*2
    print("pA: ", pA, "\tpB: ", pB, "\tN: ", n)
    return n


N_PATIENTS_DICT = {}
for pair in TRUE_P_PAIRS:
    pApB = [tok.split("=")[1] for tok in pair.split("_")]
    pA = float(pApB[0])
    pB = float(pApB[1])
    N_PATIENTS_DICT["_".join(pApB)] = compute_n_patients(pA,pB)


def get_n_patients(wc):
    return N_PATIENTS_DICT["{}_{}".format(wc["pA"], wc["pB"])]


# /NUMBER OF PATIENTS
###############################

##########################################
# SIMULATION STUDY RULES
##########################################

rule all:
    input:
        # Simulation study results
        expand(path.join(FIG_DIR, "histories", "design=blockraropt", "{probs}_fc={fc}_bc={bc}_pr={pr}_stat={stat}.png"),
               probs=TRUE_P_PAIRS,
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               pr=PRIOR_STRENGTH,
               stat=OPT_STAT),
        expand(path.join(FIG_DIR, "histories", "design=rar", "{probs}.png"),
               probs=TRUE_P_PAIRS),
        expand(path.join(FIG_DIR, "histories", "design=blockrar", "{probs}.png"),
               probs=TRUE_P_PAIRS),
        expand(path.join(FIG_DIR, "frontiers", "{probs}_bc={bc}_pr={pr}_stat={stat}.png"),
               probs=TRUE_P_PAIRS, bc=BLOCK_COST, pr=PRIOR_STRENGTH, stat=OPT_STAT),
        expand(path.join(LATEX_COMPARISON_DIR, "fc={fc}_bc={bc}_pr={pr}_stat={stat}_scenario={scen}_suffix={suff}.tex"),
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               pr=PRIOR_STRENGTH,
               stat=OPT_STAT,
               scen=["null","alternative"],
               suff=["_z",""]),
        # Trial redesign results
        expand(path.join(RED_TEST_SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_pr={pr}_stat={stat}.tsv"),
	       fc=[REDESIGN_CHOSEN_FC], bc=[REDESIGN_CHOSEN_BC], pr=[REDESIGN_PR], stat=[REDESIGN_STAT]),
        expand(path.join(RED_FIG_DIR, "swept_frontiers", "{probs}_pr={pr}_stat={stat}.png"),
	       probs=TRUE_P_PAIRS, pr=[REDESIGN_PR], stat=[REDESIGN_STAT]) 



##########################################
# SIMULATION STUDY RULES
##########################################

rule comparison_to_latex:
    input:
        src=path.join(SCRIPT_DIR, "comparison_to_latex.py"),
        comp=path.join(COMPARISON_DIR, "{params}.xlsx")
    output:
        path.join(LATEX_COMPARISON_DIR, "{params}_scenario={scen,[a-z]+}_suffix={suffix,[_a-z]*}.tex")
    shell:
        "python {input.src} {input.comp} {output} --scenarios {wildcards.scen} --suffix {wildcards.suffix}"


rule compare_designs:
    input:
        src=path.join(SCRIPT_DIR, "tabulate_comparisons.py"),
        blockraropt_tsv=path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_{bro_params}.tsv"),
        other_tsvs=expand(path.join(SCORE_DIR, "design={des}", "fc={{fc}}_bc={{bc}}_.tsv"),
                          des=BASELINES)
    output:
        path.join(COMPARISON_DIR, "fc={fc,[.0-9]+}_bc={bc,[.0-9]+}_{bro_params}.xlsx")
    shell:
        "python {input.src} --score_tsvs {input.other_tsvs} {input.blockraropt_tsv} --identifiers {BASELINES} blockraropt --output_tsv {output}" 



def get_kv_pairs(wc):
    return " ".join(wc["other_params"].strip("_").split("_"))


rule plot_fc_bc:
    input:
        src=path.join(SCRIPT_DIR, "plot_fc_bc.py"),
        tsvs=expand(path.join(SCORE_DIR, "design={{design}}","fc={fc}_bc={bc}_{{other_params}}.tsv"),
	            fc=FAILURE_COST,
		    bc=BLOCK_COST)
    output:
        path.join(FIG_DIR, "fc_bc", "design={design}", "score={score}_pA={pA}_pB={pB,[.0-9]+}_{other_params,[_.=0-9a-z]*}.png")
    shell:
        "python {input.src} {wildcards.design} {wildcards.pA} {wildcards.pB} {wildcards.score} {output} --score_tsvs {input.tsvs}"



rule plot_frontiers:
    input:
        src=path.join(SCRIPT_DIR, "plot_frontier.py"),
        bl_scores=expand(path.join(SCORE_DIR,"design={design}","fc={fc}_bc={{bc}}_.tsv"), design=BASELINES, fc=FAILURE_COST),
	bro_scores=expand(path.join(SCORE_DIR, "design=blockraropt","fc={fc}_bc={{bc}}_pr={{pr}}_stat={{stat}}.tsv"), fc=FAILURE_COST)
    output:
        os.path.join(FIG_DIR, "frontiers", "pA={pA}_pB={pB}_bc={bc}_pr={pr}_stat={stat}.png")
    shell:
        "python {input.src} {wildcards.pA} {wildcards.pB} {SCORE_DIR} {output} --bc {wildcards.bc} --pr {wildcards.pr} --stat {wildcards.stat} --baselines {BASELINES}"


rule plot_histories:
    input:
        js=path.join(HISTORY_DIR, "design={design}", "{other_params}.json"),
        src=path.join(SCRIPT_DIR, "plot_histories.py")
    output:
        path.join(FIG_DIR, "histories", "design={design}", "{other_params}.png")
    shell:
        "python {input.src} {input.js} {wildcards.design} {output}"


rule score_and_aggregate_blockraropt:
    input:
        src=path.join(SCRIPT_DIR, "score_and_aggregate.py"),
        summ=expand(path.join(SUMMARY_DIR, "design=blockraropt", "{probs}_fc={{fc}}_bc={{bc}}_pr={{pr}}_stat={{stat}}.json"), 
                    probs=TRUE_P_PAIRS)
    output:
        path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_pr={pr}_stat={stat}.tsv")
    shell:
        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"


rule score_and_aggregate_default:
    input:
        src=path.join(SCRIPT_DIR, "score_and_aggregate.py"),
        summ=expand(path.join(SUMMARY_DIR, "design={{des}}", "{probs}.json"), 
                    probs=TRUE_P_PAIRS)
    output:
        path.join(SCORE_DIR, "design={des}", "fc={fc}_bc={bc}_.tsv")
    shell:
        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"


rule summarize_rar_histories:
    input:
        src=path.join(SCRIPT_DIR, "summarize_histories.R"),
        hist=path.join(HISTORY_DIR, "design=rar", "pA={pA}_pB={pB}.json")
    output:
        path.join(SUMMARY_DIR, "design=rar", "pA={pA}_pB={pB,[.0-9]+}.json")
    shell:
        "Rscript {input.src} {input.hist} {output} --stratum_frac 0.5 --pa {wildcards.pA} --pb {wildcards.pB}"


rule summarize_histories:
    input:
        src=path.join(SCRIPT_DIR, "summarize_histories.R"),
        hist=path.join(HISTORY_DIR, "design={des}", "pA={pA}_pB={pB}_{other_params}.json")
    output:
        path.join(SUMMARY_DIR, "design={des}", "pA={pA}_pB={pB,[.0-9]+}_{other_params}.json")
    shell:
        "Rscript {input.src} {input.hist} {output} --pa {wildcards.pA} --pb {wildcards.pB}"


rule summarize_histories_default:
    input:
        src=path.join(SCRIPT_DIR, "summarize_histories.R"),
        hist=path.join(HISTORY_DIR, "design={des}", "pA={pA}_pB={pB}.json")
    output:
        path.join(SUMMARY_DIR, "design={des,(traditional|blockrar)}", "pA={pA}_pB={pB,[.0-9]+}.json")
    shell:
        "Rscript {input.src} {input.hist} {output} --pa {wildcards.pA} --pb {wildcards.pB}"


rule simulate_traditional_design:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
    output:
        path.join(HISTORY_DIR, "design=traditional", "pA={pA}_pB={pB}.json")
    params:
        n_pat=get_n_patients
    shell:
        "Rscript {input.rscript} {params.n_pat} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output}"


rule simulate_rar_design:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
    output:
        path.join(HISTORY_DIR, "design=rar", "pA={pA}_pB={pB}.json")
    params:
        n_pat=get_n_patients
    shell:
        "Rscript {input.rscript} {params.n_pat} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output} --design rar"


rule simulate_blockrar_design:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
    output:
        path.join(HISTORY_DIR, "design=blockrar", "pA={pA}_pB={pB}.json")
    params:
        n_pat=get_n_patients
    shell:
        "Rscript {input.rscript} {params.n_pat} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output} --design blockrar --n_blocks {BLOCKRAR_N_BLOCKS}"


def get_policy_db(wc):
    return path.join(POLICY_DB_DIR, "pat="+str(get_n_patients(wc))+"_{design}_stat={stat}.sqlite")

rule simulate_blockraropt_design:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
        db=get_policy_db
    output:
        path.join(HISTORY_DIR, "design=blockraropt", "pA={pA}_pB={pB,[.0-9]+}_{design}_stat={stat}.json")
    params:
        n_pat=get_n_patients
    shell:
        "Rscript {input.rscript} {params.n_pat} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output} --design blockraropt --blockraropt_db {input.db}"


def compute_block_incr(wc):
    pat = int(wc["pat"])
    inc = 2
    if pat >= 92: # 100:
        inc = 4
    elif pat > 128:
        inc = 8
    return inc 

def compute_min_blocksize(wc):
    pat = int(wc["pat"])
    return max(4, pat//8)

rule precompute_blockraropt:
    input:
        path.join(SCRIPT_DIR, "trialmdp_wrapper.R")
    output:
        path.join(POLICY_DB_DIR, "pat={pat}_fc={fc}_bc={bc}_pr={pr}_stat={stat}.sqlite")
    params:
        min_blocksize=compute_min_blocksize,
        block_incr=compute_block_incr
    shell:
        "Rscript {input} {wildcards.pat} {wildcards.fc} {wildcards.bc} {output} {params.min_blocksize} {params.block_incr} --alpha_a {wildcards.pr} --beta_a {wildcards.pr} --alpha_b {wildcards.pr} --beta_b {wildcards.pr} --test_statistic {wildcards.stat} --act_n {ACT_N}"




##########################################
# TRIAL REDESIGN RULES
##########################################

REDESIGN_BC_STR = " ".join("{}".format(bc) for bc in REDESIGN_BC)

rule redesign_test_score_and_aggregate:
    input:
        src=path.join(SCRIPT_DIR, "score_and_aggregate.py"),
        summ=expand(path.join(RED_SUMMARY_DIR, "design=blockraropt", "{probs}_fc={{fc}}_bc={{bc}}_pr={{pr}}_stat={{stat}}.json"), 
                    probs=REDESIGN_P_PAIRS)
    output:
        path.join(RED_TEST_SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_pr={pr}_stat={stat}.tsv")
    shell:
        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"

#rule redesign_chosen_simulations:
#    input:
#        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
#        db=path.join(RED_POLICY_DB_DIR, "fc={fc}_bc={bc}_pr={pr}_stat={stat}.sqlite")
#    output:
#        path.join(RED_HISTORY_DIR, "design=blockraropt", "pA={pA}_pB={pB,[.0-9]+}_fc={fc}_bc={bc}_pr={pr}_stat={stat}.json")
#    shell:
#        "Rscript {input.rscript} {REDESIGN_N} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output} --design blockraropt --blockraropt_db {input.db}"


rule redesign_sweep_plot_frontiers:
    input:
        src=path.join(SCRIPT_DIR, "plot_frontier.py"),
	bro_scores=expand(path.join(RED_SCORE_DIR, "design=blockraropt","fc={fc}_bc={bc}_pr={{pr}}_stat={{stat}}.tsv"), fc=REDESIGN_FC, bc=REDESIGN_BC)
    output:
        os.path.join(RED_FIG_DIR, "swept_frontiers", "pA={pA}_pB={pB}_pr={pr}_stat={stat}.png")
    shell:
        "python {input.src} {wildcards.pA} {wildcards.pB} {RED_SCORE_DIR} {output} --bc {REDESIGN_BC_STR} --pr {REDESIGN_PR} --stat {REDESIGN_STAT} --baselines --xl 0.0 --xu 0.6"


rule redesign_sweep_score_and_aggregate:
    input:
        src=path.join(SCRIPT_DIR, "score_and_aggregate.py"),
        summ=expand(path.join(RED_SUMMARY_DIR, "design=blockraropt", "{probs}_fc={{fc}}_bc={{bc}}_pr={{pr}}_stat={{stat}}.json"), 
                    probs=TRUE_P_PAIRS)
    output:
        path.join(RED_SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_pr={pr}_stat={stat}.tsv")
    shell:
        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"


rule redesign_summarize_histories:
    input:
        src=path.join(SCRIPT_DIR, "summarize_histories.R"),
        hist=path.join(RED_HISTORY_DIR, "design={des}", "pA={pA}_pB={pB}_{other_params}.json")
    output:
        path.join(RED_SUMMARY_DIR, "design={des}", "pA={pA}_pB={pB,[.0-9]+}_{other_params}.json")
    shell:
        "Rscript {input.src} {input.hist} {output} --pa {wildcards.pA} --pb {wildcards.pB}"


rule redesign_sweep_simulation:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
        db=path.join(RED_POLICY_DB_DIR, "fc={fc}_bc={bc}_pr={pr}_stat={stat}.sqlite")
    output:
        path.join(RED_HISTORY_DIR, "design=blockraropt", "pA={pA}_pB={pB,[.0-9]+}_fc={fc}_bc={bc}_pr={pr}_stat={stat}.json")
    shell:
        "Rscript {input.rscript} {REDESIGN_N} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output} --design blockraropt --blockraropt_db {input.db}"


rule redesign_precompute_blockraropt:
    input:
        path.join(SCRIPT_DIR, "trialmdp_wrapper.R")
    output:
        path.join(RED_POLICY_DB_DIR, "fc={fc}_bc={bc}_pr={pr}_stat={stat}.sqlite")
    shell:
        "Rscript {input} {REDESIGN_N} {wildcards.fc} {wildcards.bc} {output} {REDESIGN_MINSIZE} {REDESIGN_INCR} --alpha_a {wildcards.pr} --beta_a {wildcards.pr} --alpha_b {wildcards.pr} --beta_b {wildcards.pr} --test_statistic {wildcards.stat} --act_n {ACT_N}"




