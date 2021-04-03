

from os import path
from scipy.stats import norm

configfile: "config.yaml"

SCRIPT_DIR = "scripts"

RESULT_DIR = config["result_dir"]
POLICY_DB_DIR = path.join(RESULT_DIR, "precomputed_policies")
HISTORY_DIR = path.join(RESULT_DIR, "histories")
SUMMARY_DIR = path.join(RESULT_DIR, "summaries") 
SCORE_DIR = path.join(RESULT_DIR, "scores")
COMPARISON_DIR = path.join(RESULT_DIR, "comparisons")
BLOCKRAROPT_COMPARISON_DIR = path.join(RESULT_DIR, "blockraropt_comparisons")
FIG_DIR = path.join(RESULT_DIR, "figures")

SIM_PARAMS = config["simulation_params"]
N_SAMPLES = SIM_PARAMS["n_samples"]
TRUE_PA = SIM_PARAMS["true_pa"]
TRUE_PB = SIM_PARAMS["true_pb"]
NULL_HYP_N = SIM_PARAMS["null_hypothesis_n"]
TRUE_P_PAIRS = ["pA={}_pB={}".format(pa, pb) for pa, pb in zip(TRUE_PA,TRUE_PB)]
SIM_ALPHA = SIM_PARAMS["alpha"]
SIM_BETA = SIM_PARAMS["beta"]

###############################
# COMPUTE NUMBER OF PATIENTS
def compute_n_patients(pA, pB, alpha=SIM_ALPHA, beta=SIM_BETA):
    
    if pA == pB:
        # "Null" case
	#delta = null_delta(pB)
	#pA = pB + delta
        return NULL_HYP_N

    delta = abs(pA - pB)

    var_b = pB*(1.0-pB)
    var_a = pA*(1.0-pA)

    k = (var_b / var_a)**0.5

    n_a = (norm.ppf(1.0-alpha) + norm.ppf(1.0-beta))**2.0 * (var_a + var_b/k) / (delta*delta)
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

BASELINES = config["baseline_designs"]

BLOCKRAROPT_PARAMS = config["blockraropt_params"]
FAILURE_COST = BLOCKRAROPT_PARAMS["failure_cost"]
ACT_N = BLOCKRAROPT_PARAMS["act_n"]
BLOCK_COST = BLOCKRAROPT_PARAMS["block_cost"]
PRIOR_STRENGTH = BLOCKRAROPT_PARAMS["prior_strength"]
OPT_STAT = BLOCKRAROPT_PARAMS["stat"]

BLOCKRAR_PARAMS = config["blockrar_params"]
BLOCKRAR_N_BLOCKS = BLOCKRAR_PARAMS["n_blocks"]


rule all:
    input:
	    #        expand(path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_pr={pr}_stat={stat}.tsv"),
	    #               fc=FAILURE_COST,
	    #               bc=BLOCK_COST,
	    #               pr=PRIOR_STRENGTH,
	    #               stat=OPT_STAT
	    #               ),
	    #        expand(path.join(SCORE_DIR, "design=blockrar", "fc={fc}_bc={bc}_.tsv"),
	    #               fc=FAILURE_COST,
	    #               bc=BLOCK_COST),
	    #        expand(path.join(SCORE_DIR, "design=rar", "fc={fc}_bc={bc}_.tsv"),
	    #               fc=FAILURE_COST,
	    #               bc=BLOCK_COST,
	    #               ),
	    #        expand(path.join(SCORE_DIR, "design=traditional", "fc={fc}_bc={bc}_.tsv"),
	    #               fc=FAILURE_COST,
	    #               bc=BLOCK_COST,
	    #               ),
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
        expand(path.join(COMPARISON_DIR, "fc={fc}_bc={bc}_pr={pr}_stat={stat}.xlsx"),
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               pr=PRIOR_STRENGTH,
               stat=OPT_STAT),
	#        expand(path.join(FIG_DIR, "fc_bc", "design=blockraropt", "score={score}_{probs}_pr={pr}_stat={stat}.png"),
	#               score=["excess_failures", "blocks", "cmh_2s", "utility_cmh"],
	#               probs=TRUE_P_PAIRS,
	#               pr=PRIOR_STRENGTH,
	#               stat=OPT_STAT),
	#        expand(path.join(FIG_DIR, "fc_bc", "design=traditional", "score={score}_{probs}_.png"),
	#               score=["excess_failures", "utility_cmh"],
	#               probs=TRUE_P_PAIRS,
	#               pr=PRIOR_STRENGTH,
	#               stat=OPT_STAT)
        path.join(BLOCKRAROPT_COMPARISON_DIR, "comparisons.xlsx")


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


rule compare_blockraropt:
    input:
        src=path.join(SCRIPT_DIR, "compare_blockraropt.py"),
        tsvs=expand(path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_pr={pr}_stat={stat}.tsv"),
			      fc=FAILURE_COST, bc=BLOCK_COST,
			      pr=PRIOR_STRENGTH, stat=OPT_STAT)
    output:
        path.join(BLOCKRAROPT_COMPARISON_DIR, "comparisons.xlsx")
    shell:
        "python {input.src} --score_tsvs {input.tsvs} --output_tsv {output}" 


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
        summ=expand(path.join(SUMMARY_DIR, "design=blockraropt", "{probs}_fc={{fc}}_bc={{bc}}_stat={{stat}}.json"), 
                    probs=TRUE_P_PAIRS)
    output:
        path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_stat={stat}.tsv")
    shell:
        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"


#rule score_and_aggregate_blockrar:
#    input:
#        src=path.join(SCRIPT_DIR, "score_and_aggregate.py"),
#        summ=expand(path.join(SUMMARY_DIR, "design=blockrar", "{probs}_nblocks={{nblocks}}.json"), 
#                    probs=TRUE_P_PAIRS)
#    output:
#        path.join(SCORE_DIR, "design=blockrar", "fc={fc}_bc={bc}_nblocks={nblocks}.tsv")
#    shell:
#        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"


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
        "Rscript {input.src} {input.hist} {output} --stratum_frac 0.25"


rule summarize_histories:
    input:
        src=path.join(SCRIPT_DIR, "summarize_histories.R"),
        hist=path.join(HISTORY_DIR, "design={des}", "pA={pA}_pB={pB}_{other_params}.json")
    output:
        path.join(SUMMARY_DIR, "design={des}", "pA={pA}_pB={pB,[.0-9]+}_{other_params}.json")
    shell:
        "Rscript {input.src} {input.hist} {output}"


rule summarize_histories_default:
    input:
        src=path.join(SCRIPT_DIR, "summarize_histories.R"),
        hist=path.join(HISTORY_DIR, "design={des}", "pA={pA}_pB={pB}.json")
    output:
        path.join(SUMMARY_DIR, "design={des,(traditional|blockrar)}", "pA={pA}_pB={pB,[.0-9]+}.json")
    shell:
        "Rscript {input.src} {input.hist} {output}"


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
    if pat > 100:
        inc = 4
    elif pat > 128:
        inc = 8
    return inc 


rule precompute_blockraropt:
    input:
        path.join(SCRIPT_DIR, "blockraropt_wrapper.R")
    output:
        path.join(POLICY_DB_DIR, "pat={pat}_fc={fc}_bc={bc}_pr={pr}_stat={stat}.sqlite")
    params:
        block_incr=compute_block_incr
    shell:
        "Rscript {input} {wildcards.pat} {params.block_incr} {wildcards.fc} {wildcards.bc} {output} --alpha_a {wildcards.pr} --beta_a {wildcards.pr} --alpha_b {wildcards.pr} --beta_b {wildcards.pr} --test_statistic {wildcards.stat} --act_n {ACT_N}"


