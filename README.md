# Example snakemake workflow for performing dominating set differential abundance analysis

Dominating set differential abundance analysis allows you to perform differential abundacne analysis directly on an assembly graph.
This approach was developed in the repository [dib-lab/2020-ibd](https://github.com/dib-lab/2020-ibd) and is explained in [taylorreiter/2021-paper-ibd](https://github.com/taylorreiter/2021-paper-ibd).

This workflow takes as input raw metagenome sequences from two groups of metagenomes (e.g. case vs. control). 
It generates a taxonomic profiles for each metagenome, and selects the species that are present in some threshold of metagenomes (by default, 100%) for dominating set differential abundance analysis.

Abreviations:
+ **dda**: dominating set differential abundance
+ **CD**: Crohn's disease
+ **ccs**: corncobs; refers to the output of the differential abundance testing software corncob.

## Running the example workflow

This is a snakemake workflow that uses conda to manage software installations.

You will need `conda` to be installed in order to run this workflow.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

First, clone the repository to your machine.

```
git clone https://github.com/taylorreiter/2022-dominating-set-differential-abundance-example.git
```

`cd` into the repo and set up the conda environment in which the workflow will be run.

```
cd 2022-dominating-set-differential-abundance-example
conda env create --name dda --file environment.yml
conda activate dda
```

Obtain the test data by downloading and untarring it.
Note the test data _must_ be in the folder `inputs/mgx_raw` for the workflow to run.

```
wget -O inputs/mgx_raw.tar.gz https://osf.io/cfkjz/download
cd inputs
tar xf mgx_raw.tar.gz
```

`cd` back into the main folder of the repo and try a dry run of the workflow:
```
cd ..
snakemake -n
```

To run the workflow, you can use:
```
snakemake --use-conda --rerun-incomplete -j 1
```

If you're working on a cluster like slurm, you can use the following command to submit each rule as a job.
You'll need to update information about the partition to match your cluster.

```
snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t {resources.time_min} -J dda -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

The test data is fairly small, but the workflow still takes ~half a day to run from start to finish.


## Setting up the workflow to run with your data


## Interpretting the output files from the workflow

The final set of results is currently written to the folder `outputs/metapangenome_sgc_catlases_corncob`. 
There are two files, `*all_ccs.tsv` and `*_sig_ccs.tsv`.   
These files record all results from the differential abundance analysis (`*all_ccs.tsv`) and results that are significant after bonferonni p value correction (`*_sig_ccs.tsv`). 
Note that the test data does not produce significant results, so these files are blank other than the column names.

## Changes to the workflow you may want to make when running on your data


## Future work

This workflow is currently missing the annotation portion.
This arm of the workflow will be added soon, and it will be based upon [this snakefile](https://github.com/dib-lab/2020-ibd/blob/master/annotate_metapangenome_species_graphs.snakefile).

I plan to provide more complete documentation in the future as well.
In the meantime, many of the rules in Snakefile have docstrings that explain what the rules do, possible implementation trade offs, etc.

I also plan to add benchmarking directives to each rule to automatically capture the RAM, CPU, and run time for each step in the workflow.

@taylorreiter 06/2022
