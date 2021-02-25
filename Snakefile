

from os import path

configfile: "config.yaml"

SCRIPT_DIR = "scripts"

RESULT_DIR = config["result_dir"]
POLICY_DB_DIR = path.join(RESULT_DIR, "precomputed_policies")
HISTORY_DIR = path.join(RESULT_DIR, "histories")
SUMMARY_DIR = path.join(RESULT_DIR, "summaries") 
SCORE_DIR = path.join(RESULT_DIR, "scores")
COMPARISON_DIR = path.join(RESULT_DIR, "comparisons")
FIG_DIR = path.join(RESULT_DIR, "figures")

SIM_PARAMS = config["simulation_params"]
N_SAMPLES = SIM_PARAMS["n_samples"]
TRUE_P_RANGE = SIM_PARAMS["true_p_range"]
TRUE_P_PAIRS = ["pA={}_pB={}".format(pa, pb) for i, pa in enumerate(TRUE_P_RANGE) for pb in TRUE_P_RANGE[:(i+1)] ]

###############################
# COMPUTE NUMBER OF PATIENTS
# A stand-in until Thevaa provides some code
def compute_n_patients(pA, pB, target_power=0.8, target_t1=0.05):
    return 128

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
BLOCK_INCR = BLOCKRAROPT_PARAMS["block_increment_frac"]
FAILURE_COST = BLOCKRAROPT_PARAMS["failure_cost"]
BLOCK_COST = BLOCKRAROPT_PARAMS["block_cost"]
PRIOR_STRENGTH = BLOCKRAROPT_PARAMS["prior_strength"]
OPT_STAT = BLOCKRAROPT_PARAMS["stat"]

BLOCKRAR_PARAMS = config["blockrar_params"]
BLOCKRAR_N_BLOCKS = BLOCKRAR_PARAMS["n_blocks"]


rule all:
    input:
        expand(path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_inc={inc}_pr={pr}_stat={stat}.tsv"),
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               inc=BLOCK_INCR,
               pr=PRIOR_STRENGTH,
               stat=OPT_STAT
               ),
        expand(path.join(SCORE_DIR, "design=blockrar", "fc={fc}_bc={bc}.tsv"),
               fc=FAILURE_COST,
               bc=BLOCK_COST),
        expand(path.join(SCORE_DIR, "design=rar", "fc={fc}_bc={bc}.tsv"),
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               ),
        expand(path.join(SCORE_DIR, "design=traditional", "fc={fc}_bc={bc}.tsv"),
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               ),
        expand(path.join(FIG_DIR, "histories", "design=blockraropt", "{probs}_fc={fc}_bc={bc}_inc={inc}_pr={pr}_stat={stat}.png"),
               probs=TRUE_P_PAIRS,
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               inc=BLOCK_INCR,
               pr=PRIOR_STRENGTH,
               stat=OPT_STAT),
        expand(path.join(FIG_DIR, "histories", "design=rar", "{probs}.png"),
               probs=TRUE_P_PAIRS),
	expand(path.join(FIG_DIR, "histories", "design=blockrar", "{probs}.png"),
               probs=TRUE_P_PAIRS),
        expand(path.join(COMPARISON_DIR, "fc={fc}_bc={bc}_inc={inc}_pr={pr}_stat={stat}.xlsx"),
               fc=FAILURE_COST,
               bc=BLOCK_COST,
               inc=BLOCK_INCR,
               pr=PRIOR_STRENGTH,
               stat=OPT_STAT)


rule compare_designs:
    input:
        src=path.join(SCRIPT_DIR, "tabulate_comparisons.py"),
        blockraropt_tsv=path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_{bro_params}.tsv"),
        other_tsvs=expand(path.join(SCORE_DIR, "design={des}", "fc={{fc}}_bc={{bc}}.tsv"),
                          des=BASELINES)
    output:
        path.join(COMPARISON_DIR, "fc={fc,[.0-9]+}_bc={bc,[.0-9]+}_{bro_params}.xlsx")
    shell:
        "python {input.src} --score_tsvs {input.other_tsvs} {input.blockraropt_tsv} --identifiers {BASELINES} blockraropt --output_tsv {output}" 


def get_kv_pairs(wc):
    return " ".join(wc["other_params"].strip("_").split("_"))


rule plot_histories:
    input:
        js=path.join(HISTORY_DIR, "design={design}", "{other_params}.json"),
        src=path.join(SCRIPT_DIR, "plot_histories.py")
    output:
        path.join(FIG_DIR, "histories", "design={design}", "{other_params}.png")
    shell:
        "python {input.src} {input.js} {output}"


rule score_and_aggregate_blockraropt:
    input:
        src=path.join(SCRIPT_DIR, "score_and_aggregate.py"),
        summ=expand(path.join(SUMMARY_DIR, "design=blockraropt", "{probs}_fc={{fc}}_bc={{bc}}_inc={{inc}}_stat={{stat}}.json"), 
                    probs=TRUE_P_PAIRS)
    output:
        path.join(SCORE_DIR, "design=blockraropt", "fc={fc}_bc={bc}_inc={inc}_stat={stat}.tsv")
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
        path.join(SCORE_DIR, "design={des}", "fc={fc}_bc={bc}.tsv")
    shell:
        "python {input.src} --summaries {input.summ} --output_tsv {output} --fc {wildcards.fc} --bc {wildcards.bc}"


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
        path.join(SUMMARY_DIR, "design={des}", "pA={pA}_pB={pB,[.0-9]+}.json")
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
    return int(float(wc["inc"])*int(wc["pat"]))

rule precompute_blockraropt:
    input:
        path.join(SCRIPT_DIR, "blockraropt_wrapper.R")
    output:
        path.join(POLICY_DB_DIR, "pat={pat}_fc={fc}_bc={bc}_inc={inc}_pr={pr}_stat={stat}.sqlite")
    params:
        block_incr=compute_block_incr
    shell:
        "Rscript {input} {wildcards.pat} {params.block_incr} {wildcards.fc} {wildcards.bc} {output} --alpha_a {wildcards.pr} --beta_a {wildcards.pr} --alpha_b {wildcards.pr} --beta_b {wildcards.pr} --test_statistic {wildcards.stat}"


