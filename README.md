# TrialMDP-analyses
Analyses, experiments, and evaluations for the [TrialMDP method](https://github.com/dpmerrell/TrialMDP).

See the [TrialMDP manuscript on arXiv](https://arxiv.org/abs/2109.14642). 
The Snakemake workflow in this repository reproduces all of the results described in that paper.

## Setup

* Set up your R environment:
    - Install the [`TrialMDP`](https://github.com/dpmerrell/TrialMDP) package
    - Install the `optparse` package

* Set up your python 3.8 environment: 
    - run `pip install -r requirements.txt` to install the requisite python packages

## Reproducing the analyses from our manuscript

* Make sure you have sufficient compute resources.
    - The analyses entail hundreds of CPU-hours (can be parallelized by `--cores` option, see below)
    - The analyses require a few GB of memory and disk

* Call [Snakemake](https://snakemake.readthedocs.io/en/v5.1.4/index.html). 
    - On a workstation, execute `snakemake --cores CORES` (where `CORES` is replaced by the number of cores you want to devote)
    - On a cluster, execute `snakemake --profile PROFILE` (where `PROFILE` is replaced by some suitable [Snakemake profile](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles) )
    
