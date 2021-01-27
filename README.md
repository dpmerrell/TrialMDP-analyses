# blockRARopt-analyses
Analyses, experiments, and evaluations for the [blockRARopt method](https://github.com/dpmerrell/blockRARopt).

See our manuscript at [MANUSCRIPT LINK]().

## Setup

* Set up your R environment:
    - Install the `blockRARopt` package

* Set up your python 3.8 environment: 
    - run `pip install -r requirements.txt` to install the requisite python packages

## Reproducing the analyses from our manuscript

* Make sure you have sufficient compute resources.
    - The analyses entail thousands of CPU-hours
    - The analyses require several GB of memory and disk

* Call [Snakemake](https://snakemake.readthedocs.io/en/v5.1.4/index.html). 
    - On a workstation, execute `snakemake --cores CORES` (where `CORES` is replaced by the number of cores you want to devote)
    - On a cluster, execute `snakemake --profile PROFILE` (where `PROFILE` is replaced by some suitable [Snakemake profile](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles) )
    
