

from os import path

configfile: "config.yaml"

SCRIPT_DIR = "scripts"

RESULT_DIR = config["result_dir"]
POLICY_DB_DIR = path.join(RESULT_DIR, "precomputed_policies")
SIM_OUTPUT_DIR = path.join(RESULT_DIR, "sim_output")
FIG_DIR = path.join(RESULT_DIR, "figures")
SCORE_DIR = path.join(RESULT_DIR, "scores") 
TABLE_DIR = path.join(RESULT_DIR, "tabulated")

SIM_PARAMS = config["simulation_params"]
N_PATIENTS = SIM_PARAMS["n_patients"]
N_SAMPLES = SIM_PARAMS["n_samples"]
TRUE_P_RANGE = SIM_PARAMS["true_p_range"]
TRUE_P_PAIRS = ["pA={}_pB={}".format(pa, pb) for i, pa in enumerate(TRUE_P_RANGE) for pb in TRUE_P_RANGE[:(i+1)] ]

BLOCKRAROPT_PARAMS = config["blockraropt_params"]
BLOCK_INCR = BLOCKRAROPT_PARAMS["block_increment_frac"]
FAILURE_COST = BLOCKRAROPT_PARAMS["failure_cost"]
BLOCK_COST = BLOCKRAROPT_PARAMS["block_cost"]
PRIOR_STRENGTH = BLOCKRAROPT_PARAMS["prior_strength"]


rule all:
    input:
        #path.join(TABLE_DIR, "blockraropt.tsv"),
	#path.join(TABLE_DIR, "traditional.tsv"),
        expand(path.join(FIG_DIR,"heatmaps","score={score}__design=traditional_pat={pat}.png"),
               score=["prob_reject", "excess_failure_frac"],
               pat=N_PATIENTS),
        expand(path.join(FIG_DIR,"heatmaps", "score={score}__design=blockraropt_pat={pat}_fc={fc}_bc={bc}_pr={pr}.png"),
               score=["prob_reject", "excess_failure_frac"],
               pat=N_PATIENTS, fc=FAILURE_COST, bc=BLOCK_COST, pr=PRIOR_STRENGTH)
               

def get_kv_pairs(wc):
    return " ".join(wc["other_params"].strip("_").split("_"))

rule plot_probspace_heatmap:
    input:
        path.join(TABLE_DIR, "{design}.tsv")
    output:
        path.join(FIG_DIR,"heatmaps","score={score}__design={design,[a-z]+}_{other_params}.png")
    params:
        kvp=get_kv_pairs
    shell:
        "python scripts/plot_heatmap.py {input} {wildcards.score} {output} --kv_pairs design={wildcards.design} {params.kvp}"


rule tabulate_blockraropt_scores:
    input:
        jsons=expand( path.join(SCORE_DIR, "design=blockraropt", "{probs}_pat={pat}_fc={fc}_bc={bc}_inc={inc}_pr={pr}.json"), probs=TRUE_P_PAIRS, pat=N_PATIENTS, inc=BLOCK_INCR, fc=FAILURE_COST, bc=BLOCK_COST, pr=PRIOR_STRENGTH),
        src=path.join(SCRIPT_DIR, "tabulate_scores.py")
    output:
        path.join(TABLE_DIR, "blockraropt.tsv")
    shell:
        "python {input.src} {input.jsons} {output}"


rule tabulate_traditional_scores:
    input:
        jsons=expand( path.join(SCORE_DIR, "design=traditional", "{probs}_pat={pat}_fc={fc}_bc={bc}.json"), probs=TRUE_P_PAIRS, pat=N_PATIENTS, fc=FAILURE_COST, bc=BLOCK_COST),
        src=path.join(SCRIPT_DIR, "tabulate_scores.py")
    output:
        path.join(TABLE_DIR, "traditional.tsv")
    shell:
        "python {input.src} {input.jsons} {output}"


rule score_sim_results:
    input:
        json=path.join(SIM_OUTPUT_DIR, "design={des}", "pA={pA}_pB={pB}_pat={pat}_fc={fc}_bc={bc}{other_params}.json"),
	src=path.join(SCRIPT_DIR, "score_sim_results.py")
    output:
        path.join(SCORE_DIR, "design={des,[a-z]+}", "pA={pA}_pB={pB}_pat={pat}_fc={fc}_bc={bc,[.0-9]+}{other_params}.json")
    shell:
        "python {input.src} {input.json} {output} --true_p_a {wildcards.pA} --true_p_b {wildcards.pB} --failure_cost {wildcards.fc} --block_cost {wildcards.bc}"


rule simulate_traditional_design:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
    output:
        path.join(SIM_OUTPUT_DIR, "design=traditional", "pA={pA}_pB={pB}_pat={pat,[0-9]+}{other_params}.json")
    shell:
        "Rscript {input.rscript} {wildcards.pat} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output}"


rule simulate_blockraropt_design:
    input:
        rscript=path.join(SCRIPT_DIR, "simulation_wrapper.R"),
        db=path.join(POLICY_DB_DIR, "pat={pat}_{design}.sqlite")
    output:
        path.join(SIM_OUTPUT_DIR, "design=blockraropt", "pA={pA}_pB={pB}_pat={pat,[0-9]+}_{design}.json")
    shell:
        "Rscript {input.rscript} {wildcards.pat} {N_SAMPLES} {wildcards.pA} {wildcards.pB} {output} --design blockraropt --blockraropt_db {input.db}"

def compute_block_incr(wc):
    return int(float(wc["inc"])*int(wc["pat"]))

rule precompute_blockraropt:
    input:
        path.join(SCRIPT_DIR, "blockraropt_wrapper.R")
    output:
        path.join(POLICY_DB_DIR, "pat={pat}_fc={fc}_bc={bc}_inc={inc}_pr={pr}.sqlite")
    params:
        block_incr=compute_block_incr
    shell:
        "Rscript {input} {wildcards.pat} {params.block_incr} {wildcards.fc} {wildcards.bc} {output} --alpha_a {wildcards.pr} --beta_a {wildcards.pr} --alpha_b {wildcards.pr} --beta_b {wildcards.pr}"


