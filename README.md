# Example snakemake workflow for performing dominating set differential abundance analysis

Dominating set differential abundance analysis allows you to perform differential abundacne analysis directly on an assembly graph.
This approach was developed in the repository [dib-lab/2020-ibd](https://github.com/dib-lab/2020-ibd) and is explained in [taylorreiter/2021-paper-ibd](https://github.com/taylorreiter/2021-paper-ibd).

This workflow takes as input raw metagenome sequences from two groups of metagenomes (e.g. case vs. control). 
It generates a taxonomic profiles for each metagenome, and selects the species that are present in some threshold of metagenomes (by default, 100%) for dominating set differential abundance analysis.

## Running the example workflow


## Setting up the workflow to run with your data


## Changes to the workflow you may want to make when running on your data


## Future work

This workflow is currently missing the annotation portion.
This arm of the workflow will be added soon, and it will be based upon [this snakefile](https://github.com/dib-lab/2020-ibd/blob/master/annotate_metapangenome_species_graphs.snakefile).
