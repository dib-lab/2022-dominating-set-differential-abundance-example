# Example snakemake workflow for performing dominating set differential abundance analysis

Dominating set differential abundance analysis allows you to perform differential abundance analysis directly on an assembly graph.
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

Obtain the test data by downloading and decompressing it.
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


## Interpretting the output files from the workflow

The final set of results is currently written to the folder `outputs/metapangenome_sgc_catlases_corncob`. 
There are two files, `*all_ccs.tsv` and `*_sig_ccs.tsv`.   
These files record all results from the differential abundance analysis (`*all_ccs.tsv`) and results that are significant after Bonferroni p value correction (`*_sig_ccs.tsv`). 
Note that the test data does not produce significant results, so these files are blank other than the column names.

The top of `*_all_ccs.tsv` looks like this:

```
mu      estimate        standard_error  t_value p_value domset_piece    separation_in_abund_model
mu.(Intercept)  -16.774977027348083     1.1383586548403701      -14.736108831799067     2.6013395959754513e-5   1   none
mu.varCD        -0.6215969647914198     1.4488716929991103      -0.4290214018225018     0.6857641221039749      1   none
phi.(Intercept) -15.169795834455815     1.5083054795318949      -10.057508933245927     1.6630950764093052e-4   1   none
mu.(Intercept)  -9.22997784200045       0.6893660903158337      -13.38908015880463      4.158510923125367e-5    2   none
mu.varCD        0.5406446864591221      0.6605571110385722      0.818467740978648       0.45032188152927166     2   none
phi.(Intercept) -8.044757987704175      0.689533175199002       -11.666962920794095     8.126876918271136e-5    2   none
mu.(Intercept)  -7.351986768920887      0.43085265156051716     -17.063807643500677     1.2649542116123428e-5   3   none
mu.varCD        -1.1623261554542317     0.6698625144098324      -1.7351712186466706     0.14323359629040192     3   none
phi.(Intercept) -7.526117770129453      0.6342089717889471      -11.866936774640905     7.48382217116417e-5     3   none
```

## Changes to the workflow you may want to make when running on your data

This workflow has _not_ yet been adapted to make re-running it on new data easy, but it is doable with a few changes.

1. **Updating the workflow to run on your own data**: To run the workflow on your own data, you'll need to make a few changes. 
    
    a. First, you'll need a new `inputs/metadata.csv` file. This file looks like this:
    ```
    sample,var
    PSM7J199,CD
    PSM7J1BJ,CD
    MSM9VZMM,nonIBD
    MSM9VZL5,nonIBD
    p8775mo1,CD
    p8816mo1,CD
    p9220mo1,nonIBD
    p9223mo1,nonIBD
    ```
    It's a csv file with the columns `sample` and `var`.
    The `sample` column records the basename of your sample, minus the suffix `_R*.fq.gz`.
    The `var` column records the two-level variable that will be compared between samples. 
    
    b. Next, you'll need to place your own paired-end fastq files in the `inputs/mgx_raw` folder. 
    The files must be named with the same sample basename as is indicated in `inputs/metadata.csv` and must end with `_R*.fq.gz`.
    
    c. Lastly, line 54 of `scripts/corncob_dda.R` needs to be updated.
    This line encodes the levels of the `var` that will be compared, and these should exactly match the two names in `inputs/metadata.csv` in the `var` column.
    The base `var` (e.g. controls) should come first.
2. **Updating the workflow to run a more complex model**: The only files that needs to be changed to enable a more complex model design are `scripts/corncob_dda.R` and `inputs/metadata.csv`. The metadata file should include all variables that you would like to include in your model, with each encoded as a new column. The `scripts/corncob_dda.R` file will need more substantial changes. We plan to add more documentation on this later, but in the meantime for inspiration, see [here](https://github.com/dib-lab/2020-ibd/blob/master/scripts/corncob_dda.R).

Note we have not tested these suggested changes for adapting to new data...buyer beware! 
We think this is everything that is needed, but may have missed something.

## Future work

This workflow is currently missing the annotation portion.
This arm of the workflow will be added soon, and it will be based upon [this snakefile](https://github.com/dib-lab/2020-ibd/blob/master/annotate_metapangenome_species_graphs.snakefile).

I plan to provide more complete documentation in the future as well.
In the meantime, many of the rules in Snakefile have docstrings that explain what the rules do, possible implementation trade offs, etc.

I also plan to add benchmarking directives to each rule to automatically capture the RAM, CPU, and run time for each step in the workflow.

@taylorreiter 06/2022
