# Example snakemake workflow for performing dominating set differential abundance analysis

Dominating set differential abundance analysis allows you to perform differential abundance analysis directly on an assembly graph.
This approach was developed in the repository [dib-lab/2020-ibd](https://github.com/dib-lab/2020-ibd) and is explained in [dib-lab/2021-paper-ibd](https://github.com/dib-lab/2021-paper-ibd).

As input, this workflow takes raw metagenome sequences from two groups of metagenomes (e.g. case vs. control). 
It generates a taxonomic profile for each metagenome, and selects the species that are present in some threshold of metagenomes (by default, 100%) for dominating set differential abundance analysis.

Abbreviations:
+ **dda**: dominating set differential abundance
+ **CD**: Crohn's disease
+ **ccs**: corncobs; refers to the output of the differential abundance testing software [corncob](https://github.com/bryandmartin/corncob).
+ **GTDB**: [Genome taxonomy database](https://gtdb.ecogenomic.org/)

## Running the example workflow

### Installation & setup
This is a snakemake workflow that uses conda to manage software installations.

You will need `conda` to be installed in order to run this workflow.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

After installing `conda`, clone the repository to your machine.

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

`cd` back into the main folder of the repo, where we will run the actual workflow.

```
cd ..
```

### Running the workflow

The workflow is written in three parts, `00_select_query_species_for_dda.snakefile`, `01_perform_dda.snakefile` and `02_annotate_dda.snakefile`. 
`00_select_query_species_for_dda.snakefile` needs to be run and completed before `01_perform_dda.snakefile`, and  `01_perform_dda.snakefile` needs to be completed before `02_annotate_dda.snakefile` is ran.

To run the workflow, you can use:
```
snakemake -s 00_select_query_species_for_dda.snakefile --use-conda --rerun-incomplete -j 1
```

When this snakefile has finished running, you can then run the next snakefile:
```
snakemake -s 01_perform_dda.snakefile --use-conda --rerun-incomplete -j 1
```

And then the final snakefile:
```
snakemake -s 02_annotate_dda.snakefile --use-conda --rerun-incomplete -j 1
```

The test data is fairly small, but the workflow still takes ~half a day to run from start to finish.

### Running the workflow on a cluster

If you're working on a cluster/remote computer with a job scheduler (e.g. slurm), you can use the following command to submit each rule as a job.
You'll need to update information about the partition to match your cluster.

```
snakemake -s 00_select_query_species_for_dda.snakefile -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t {resources.time_min} -J dda -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

### Dry runs

A dry run let's you see what rules still need to be run, and how many times each will be run.
Each dry run will only work after the preceding snakefile is finished running as there are file dependencies between the snakefiles.

You can start by trying a dry run of the first part of the workflow:

```
snakemake -s 00_select_query_species_for_dda.snakefile -n
```

The bottom of the dry run should produce an output that looks like the one shown below.
It details all of the rules that need to be run, the number of times each rule will be run, and the threads that will be allocated to each run.
The rules are listed in alphabetical order, not in order of execution.

```
Building DAG of jobs...
Job stats:
job                                               count    min threads    max threads
----------------------------------------------  -------  -------------  -------------
all                                                   1              1              1
download_human_host_db                                1              1              1
download_sourmash_gather_database                     1              1              1
download_sourmash_gather_database_lineages            1              1              1
gunzip_sourmash_gather_database_lineages              1              1              1
mgx_fastp_reads                                       8              1              1
mgx_kmer_trim_reads                                   8              1              1
mgx_remove_host_reads                                 8              1              1
mgx_select_query_genomes_shared_across_samples        1              1              1
mgx_sourmash_gather                                   8              1              1
mgx_sourmash_sketch                                   8              1              1
total                                                46              1              1
```

The bottom of the dry run for `01_perform_dda.snakefile` looks like:
```
Building DAG of jobs...
Job stats:
job                                                    count    min threads    max threads
---------------------------------------------------  -------  -------------  -------------
all                                                        1              1              1
corncob_for_dominating_set_differential_abund              3              1              1
diginorm_spacegraphcats_query_genomes                      3              1              1
download_query_genome                                      3              1              1
format_spacegraphcats_pangenome_catlas_abundances          3              1              1
generate_charcoal_genome_list                              1              1              1
hardtrim_spacegraphcats_query_genomes                      3              1              1
make_query_genome_info_csv                                 3              1              1
make_sgc_metapangenome_conf_files                          3              1              1
metapangeome_spacegraphcats_build                          3              1              1
mgx_make_sgc_conf                                          8              1              1
mgx_spacegraphcats_query_genomes_extract_reads             8              1              1
query_genomes_charcoal_decontaminate                       1              1              1
spacegraphcats_pangenome_catlas_estimate_abundances        3              1              1
tmp_cp_sgc_nbhds_w_lib_prefix                             24              1              1
total                                                     70              1              1
```

The bottom of the dry run for `02_annotate_dda.snakefile` looks like:
```
Job stats:
job                                                                           count    min threads    max threads
--------------------------------------------------------------------------  -------  -------------  -------------
all                                                                               1              1              1
bakta_annotate_gtdb_genomes                                                     182              1              1
bakta_download_db                                                                 1              1              1
cat_annotated_sequences                                                           3              1              1
cluster_annotated_sequences                                                       3              1              1
eggnog_annotate_clustered_sequences                                               3              1              1
eggnog_download_db                                                                1              1              1
grab_annotation_genome_accessions                                                 3              1              1
join_annotations_to_dominating_set_differential_abundance_analysis_results        3              1              1
make_metapangenome_sgc_multifasta_conf_files                                      3              1              1
sketch_metapangenome_reference                                                    3              1              1
spacegraphcats_metapangenome_catlas_cdbg_to_pieces_map                            3              1              1
spacegraphcats_metapangenome_catlas_multifasta_annotate                           3              1              1
translate_clustered_sequences_for_annotations                                     3              1              1
total                                                                           215              1              1
```

## Overview of the steps in the workflow

The workflow is written in three parts.
The first snakefile selects genomes with which to perform differential abundance analysis, the second snakefile performs the differential abundance analysis, and the third snakefile annotates the results with gene and ortholog names.
Below, we provide a more detailed overview of the steps encoded by each workflow.

`00_select_query_species_for_dda.snakefile`:
+ preprocesses the raw metagenome sequencing reads.
+ determines the taxonomic composition of each metagenome.
+ selects species on which to perform differential abundance analysis using species that are detected in 100% of samples (threshold can be decreased by the user). For each selected species, the GTDB representative genome becomes the _query genome_.

`01_perform_dda.snakefile`:
+ recovers all of the sequencing reads/variation (neighborhoods) attributable to the query genomes within each metagenome.
+ combines the genome query neighborhoods to build a _metapangenome graph_ that represents all of the variation observed across all metagenomes for a given species.
+ performs dominating set differential abundance analysis on the metapangenome graph to identify the sequences (dominating set pieces) that are more or less abundant in one group of metagenomes when compared to another.

`02_annotate_dda.snakefile`:
+ annotates the metapangenome graph using matches to genes in GTDB (rs207) genomes of the same species as the query genome.
+ joins annotations to the differentially abundant sequences.

## Interpreting the output files from the workflow

The final set of results is currently written to the folder `outputs/metapangenome_sgc_catlases_corncob`. 
There are two files, `*all_ccs.tsv` and `*_sig_ccs.tsv`.   
The first file records all results from the differential abundance analysis (`*all_ccs.tsv`) and while the second records results that are significant after Bonferroni p value correction (`*_sig_ccs.tsv`). 
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
2. **Updating the workflow to run a more complex model**: The files that need to be changed to enable a more complex model design are `scripts/corncob_dda.R` and `inputs/metadata.csv`. The metadata file should include all variables that you would like to include in your model, with each encoded as a new column. The `scripts/corncob_dda.R` file will need more substantial changes. We plan to add more documentation on this later, but in the meantime for inspiration, see [here](https://github.com/dib-lab/2020-ibd/blob/master/scripts/corncob_dda.R).

Note, we have not tested these suggested changes for adapting to new data...buyer beware! 
We think this is everything that is needed, but may have missed something.
Please file an issue if you find something to be amiss.

## Future work

I plan to provide more complete documentation in the future.
In the meantime, many of the rules in snakefiles have docstrings that explain what the rules do, possible implementation trade offs, etc.

I also plan to add benchmarking directives to each rule to automatically capture the RAM, CPU, and run time for each step in the workflow.

## Common problems when running the workflow

Most issues that pop up when running this workflow are either related to snakemake or conda.
We've documented a few common ones below -- thanks to @ccbaumler for alpha testing the workflow and documenting these problems!

### Issues related to conda software management

1. **Conda installation is not configured to use strict channel priorities:**
```
$ snakemake --use-conda --rerun-incomplete -j 1

Building DAG of jobs...
CreateCondaEnvironmentException:
Your conda installation is not configured to use strict channel priorities. This is however crucial for having robust and correct environments (for details, see https://conda-forge.org/docs/user/tipsandtricks.html). Please configure strict priorities by executing 'conda config --set channel_priority strict'.
```
This error is fixed by configuring strict priorities in conda:
```
conda config --set channel_priority strict
```

### Issues related to snakemake execution of the workflow

1. **Missing files after 5 seconds.** 
```
OSError: Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:
outputs/mgx_sgc_genome_queries/p8775mo1_k31_r1_search_oh0/GCF_008121495.1_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz
```
This problem occurs when a rule successfully finishes running, but because of file system latency, the file isn't visible to snakemake within the time period it expects to see it.
This can occur with any file output by the workflow, not just the file shown above.
The error can essentially be ignored, and the workflow restarted.
To prevent the error in the future, you can use the flag `--latency-wait` to increase the wait time. 
I've had some success with `--latency-wait 15` (15 seconds), but even that isn't always sufficient and you may need to restart the workflow even with this parameter included.

@taylorreiter 06/2022
