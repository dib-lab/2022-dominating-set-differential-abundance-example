import pandas as pd
import csv
import glob
import os

m = pd.read_csv("inputs/metadata.csv", header = 0)
SAMPLES = m['sample'].unique().tolist()

rule all:
    input: "outputs/query_genomes_from_sourmash_gather/query_genomes.csv",

########################################
## PREPROCESSING
########################################

rule mgx_fastp_reads:
    """
    fastp removes adapter sequences, removes reads that are shorter than 31 base pairs (k-mer size to build cDBG and for taxonomic discovery), and lightly quality controls reads (Phred score < 4).
    """
    input: 
        r1 = "inputs/mgx_raw/{sample}_R1.fq.gz",
        r2 = "inputs/mgx_raw/{sample}_R2.fq.gz"
    output: 
        r1 = 'outputs/mgx_fastp/{sample}_R1.fastp.fq.gz',
        r2 = 'outputs/mgx_fastp/{sample}_R2.fastp.fq.gz',
        html = 'outputs/mgx_fastp/{sample}.html',
        json = 'outputs/mgx_fastp/{sample}.json'
    conda: 'envs/fastp.yml'
    threads: 1
    resources: 
        mem_mb = 8000,
        time_min = 720
    shell:'''
    fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -q 4 --html {output.html} -j {output.json} -l 31 -c
    '''

rule download_human_host_db:
    output: "inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"
    threads: 1
    resources: 
        mem_mb = 1000,
        time_min = 30
    shell:'''
    wget -O {output} https://osf.io/84d59/download
    '''

rule mgx_remove_host_reads:
    """
    Removes host genome sequence from microbiome.
    This rule assumes a human host for a metagenome.
    If your metagenome did not come from a human but did come from a host-associated microbiome, it's probably still a good idea to try and subtract the host genome.
    While here we used a masked version of the human genome, using the genome itself as the host will still work (no promises around mitochondira/chlorophyll).
    Replace the path to the host file with the path to your host genome file.
    We've provided download scripts for a variety of common mammalian hosts in host_sequences.snakefile (mostly because we had then laying around); any genome sequence could be used here however.
    These sequences can be used as drop-in replacements for the "host" sequence indicated below.
    If your metagenome did not come from a host-associated microbiome, it's not a bad idea to keep this rule as is and remove human, as human DNA could have snuck in during DNA extraction/library prep/sequencing. 
    This method was introduced in this seqanswers post: http://seqanswers.com/forums/archive/index.php/t-42552.html
    The database was originally downloaded from this url: https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing 
    """
    input: 
        r1 = 'outputs/mgx_fastp/{sample}_R1.fastp.fq.gz',
        r2 = 'outputs/mgx_fastp/{sample}_R2.fastp.fq.gz',
        host='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    output:
        r1 = 'outputs/mgx_bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/mgx_bbduk/{sample}_R2.nohost.fq.gz',
        host_r1='outputs/mgx_bbduk/{sample}_R1.host.fq.gz',
        host_r2='outputs/mgx_bbduk/{sample}_R2.host.fq.gz'
    conda: 'envs/bbmap.yml'
    threads: 1
    resources: 
        mem_mb = 64000,
        time_min = 120
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.host}
    '''

rule mgx_kmer_trim_reads:
    """
    Trims k-mers that are probably errors. 
    The -V option (--variable-coverage) prevents elimination of low-abundance reads by only trimming low-abundance k-mers from high-abundance reads.
    K-mer trimming is an important pre-processing step before building the cDBG, as erroneous k-mers will increase the complexity of the graph.
    """
    input: 
        ancient('outputs/mgx_bbduk/{sample}_R1.nohost.fq.gz'),
        ancient('outputs/mgx_bbduk/{sample}_R2.nohost.fq.gz')
    output: "outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'envs/spacegraphcats.yml'
    threads: 1
    resources: 
        mem_mb = 64000,
        time_min = 720
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule mgx_ntcard_count_kmers_per_sample:
    input: ancient("outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz")
    output: 
        fstat = "outputs/mgx_ntcard/{sample}.fstat",
        freq = "outputs/mgx_ntcard/{sample}.freq"
    conda: 'envs/ntcard.yml'
    threads: 4
    resources:
        mem_mb=4000, 
        time_min=60
    shell:'''
    ntcard -k31 -c2000 -t {threads} -o {output.freq} {input} &> {output.fstat}
    '''

rule mgx_format_ntcard_kmer_count:
    input: fstat = expand("outputs/mgx_ntcard/{sample}.fstat", sample = SAMPLES)
    output: tsv = 'outputs/mgx_ntcard/all_kmer_count.tsv'
    conda: "envs/tidyverse.yml"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=120
    script: "scripts/format_ntcard_kmer_count.R"

##########################################################
## Determine taxonomic profile of metagenomes &
## identified species that are present in all metagenomes
##########################################################

rule mgx_sourmash_sketch:
    """
    Create a FracMinHash sketch of the quality-controlled reads.
    """
    input: ancient("outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz")
    output: "outputs/mgx_sourmash_sigs/{sample}.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        time_min=1200,
    conda: "envs/sourmash.yml"
    shell:"""
    sourmash sketch dna -p k=21,k=31,k=51,scaled=2000,abund -o {output} --name {wildcards.sample} {input}
    """

rule download_sourmash_gather_database:
    output: "inputs/gtdb-rs207.genomic-reps.dna.k31.zip"
    threads: 1
    resources:
        mem_mb = 800,
        time_min = 30
    shell:'''
    wget -O {output} https://osf.io/3a6gn/download
    '''

rule download_sourmash_gather_database_lineages:
    output: "inputs/gtdb-rs207.taxonomy.csv.gz"
    threads: 1
    resources:
        mem_mb = 800,
        time_min = 30
    shell:'''
    wget -O {output} https://osf.io/v3zmg/download
    '''

rule gunzip_sourmash_gather_database_lineages:
    input: "inputs/gtdb-rs207.taxonomy.csv.gz"
    output: "inputs/gtdb-rs207.taxonomy.csv"
    threads: 1
    resources:
        mem_mb = 800,
        time_min = 30
    shell:'''
    gunzip -c {input} > {output}
    '''

rule mgx_sourmash_gather:
    """
    Determine the taxonomic profile of each metagenome using sourmash gather.
    Compares against the GTDB rs207 reps database.
    If the organisms in your metagenome are well-represented by GTDB, the representatives database is a good choice for this step as it is fast and good enough to generate candidate query genomes for spacegraphcats and dominating set differential abundance analysis.
    If this step doesn't produce many matches, the database can be substituted for GenBank databases (built March 2022, see: https://sourmash.readthedocs.io/en/latest/databases.html).
    """
    input:
        sig=ancient("outputs/mgx_sourmash_sigs/{sample}.sig"),
        db="inputs/gtdb-rs207.genomic-reps.dna.k31.zip",
    output: 
        csv="outputs/mgx_sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv",
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 12000,
        time_min = 720
    threads: 1
    benchmark: "benchmarks/mgx/{sample}_gather.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

rule mgx_select_query_genomes_shared_across_samples:
    """
    Read in the gather results for all samples, and select query genomes that are present in some minimum fraction of samples.
    The fraction of samples a genome must be detected in is set in params.min_sample_frac, and by default is set to 0.6 (it should probably be set to 1, but these samples are so dramatically reduced in size that 0.6 made more sense for the test data set).
    We intentionally use the reps database so that there will be more shared genomes detected between samples, at the expense of the more of the sample being identifiable.
    (E.g. if E. coli is present in all samples, we'll get the same E. coli reps match instead of the best matching genome, which may be different between samples.)
    """
    input:
        gather=expand("outputs/mgx_sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv", sample = SAMPLES),
        lineages="inputs/gtdb-rs207.taxonomy.csv"
    output:
        query_genomes="outputs/query_genomes_from_sourmash_gather/query_genomes.csv",
    params: min_sample_frac = 1
    conda: 'envs/tidyverse.yml'
    resources:
        mem_mb = 4000,
        time_min = 60
    script: "scripts/mgx_select_query_genomes_shared_across_samples.R"
