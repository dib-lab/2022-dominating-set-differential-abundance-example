import pandas as pd
import csv
import glob
import os

m = pd.read_csv("inputs/metadata.csv", header = 0)
SAMPLES = m['sample'].unique().tolist()

query_genomes = pd.read_csv("outputs/query_genomes_from_sourmash_gather/query_genomes.csv", header = 0)
ACC = query_genomes['accession'].tolist()

rule all:
    input: expand("outputs/metapangenome_sgc_catlases_corncob/{acc}_sig_ccs.tsv", acc = ACC)
        

################################################################
## PREP QUERY GENOMES AND QUERY ASSEMBLY GRAPHS
################################################################

rule make_query_genome_info_csv:
    output: csvfile = 'outputs/query_genomes/{acc}.info.csv'
    conda: "envs/lxml.yml"
    resources:
        mem_mb = 8000,
        time_min = 5
    threads: 1
    shell: """
    python scripts/genbank_genomes.py {wildcards.acc} --output {output.csvfile}
    """

rule download_query_genome:
    input: csvfile = ancient('outputs/query_genomes/{acc}.info.csv')
    output: genome = "outputs/query_genomes/{acc}_genomic.fna.gz"
    resources:
        mem_mb = 500,
        time_min = 30
    threads: 1
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading genome for acc {acc}/{name} from NCBI...", file=sys.stderr)
            with open(output.genome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.genome}",
                        file=sys.stderr)


rule generate_charcoal_genome_list:
    input:  ancient(Checkpoint_GrabAccessions("outputs/query_genomes/{acc}_genomic.fna.gz"))
    output: "outputs/charcoal_conf/charcoal.genome-list.txt"
    threads: 1
    resources:
        mem_mb=500,
        time_min = 10
    shell:'''
    ls outputs/query_genomes/*gz | xargs -n 1 basename > {output} 
    '''

rule query_genomes_charcoal_decontaminate:
    """
    Because a relatively low containment index (0.01) between a query and a metagenome is necessary to recover ~20-40% of an organisms genome from the metagenome, contamination could potentially lead to recovery of many off target sequences. 
    Charcoal decontaminates the query genomes prior to running the assembly graph queries.
    """
    input:
        genomes = ancient(expand("outputs/query_genomes/{acc}_genomic.fna.gz", acc = ACC)), 
        genome_list = "outputs/charcoal_conf/charcoal.genome-list.txt",
        conf = "inputs/charcoal-conf.yml",
        genome_lineages="outputs/query_genomes_from_sourmash_gather/query_genomes.csv",
        db="inputs/gtdb-rs207.genomic-reps.dna.k31.zip",
        db_lineages="inputs/gtdb-rs207.taxonomy.csv"
    output: 
        expand("outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz", acc = ACC),
        "outputs/query_genomes_charcoal/stage1_hitlist.csv"
    resources:
        mem_mb = 64000,
        time_min = 1440
    threads: 8
    conda: "envs/charcoal.yml"
    shell:'''
    python -m charcoal run {input.conf} -j {threads} clean --nolock --latency-wait 15 --rerun-incomplete
    '''

rule mgx_make_sgc_conf:
    """
    Make configuration file that will be used for mgx_spacegraphcats_built_catlas and mgx_spacegraphcats_query_genomes_extract_reads
    """
    input: queries = expand("outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz", acc = ACC)
    output: conf = "outputs/sgc_conf/{sample}_k31_r1_conf.yml"
    resources:
        mem_mb = 500,
        time_min = 20
    threads: 1
    run:
        query_list = "\n- ".join(input.queries)
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.sample}
input_sequences:
- outputs/mgx_abundtrim/{wildcards.sample}.abundtrim.fq.gz
ksize: 31
radius: 1
paired_reads: true
search:
- {query_list}
""", file=fp)

rule mgx_spacegraphcats_query_genomes_extract_reads:
    """
    The RAM necessary to run this rule is mostly tied to the construction of the cDBG. 
    This becomes more RAM-intensive when there are more k-mers in a sample.
    As such, the RAM necessary to run this step increases with the complexity/taxonomic diversity of the sample.
    We include some suggestions based on using this approach with metagenomes from a variety of environments.
    RAM suggestions (that might be very wrong, but are probably good starting points):
    infant gut microbiome (2-20 species): 24 Gb
    mammalian gut microbiome (100-1000 species): 32-64 Gb
    ruminant gut microbiome: 500 Gb
    ocean microbiome: 1Tb
    soil microbiome: ?
    """
    input: 
        conf = ancient("outputs/sgc_conf/{sample}_k31_r1_conf.yml"),
        reads = "outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz",
        queries = expand("outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz", acc = ACC)
    output: expand("outputs/mgx_sgc_genome_queries/{{sample}}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz", acc = ACC)
    params: outdir = "outputs/mgx_sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 32000,
        time_min = 1440
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

##########################################################
## Build metapangenome graphs
##########################################################

rule diginorm_spacegraphcats_query_genomes:
    input:
        reads = expand("outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/{{acc}}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz", sample = SAMPLES)
    output: "outputs/mgx_sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    resources:
        mem_mb = 164000,
        time_min = 1440
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    cat {input.reads} | zcat | normalize-by-median.py -k 20 -C 20 -M 164e9 --gzip -o {output} -
    '''

rule hardtrim_spacegraphcats_query_genomes:
    input: "outputs/mgx_sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    output: "outputs/mgx_sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz"
    resources:
        mem_mb = 24000,
        time_min = 1440
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    trim-low-abund.py -C 4 -M 20e9 -k 31 {input} --gzip -o {output}
    '''

# POTENTIAL TODO: could add multifasta queries in here, but that would mean all of the multifasta steps would need to be completed before the metaspecies catlas could be built...
# so probably do two separate conf files. Check in testing if I need to do snakemake touch before running multifasta for it to not trigger rebuild of the catlas if using a new conf file.
rule make_sgc_metapangenome_conf_files:
    input: reads = "outputs/mgx_sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
    output: conf = "outputs/sgc_conf/{acc}_r10_conf.yml"
    resources:
        mem_mb = 500,
        time_min = 5
    threads: 1
    run:
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.acc}
input_sequences:
- {input.reads}
radius: 10
paired_reads: true
""", file=fp)

rule metapangeome_spacegraphcats_build:
    input: conf = "outputs/sgc_conf/{acc}_r10_conf.yml"
    output: 
        "outputs/metapangenome_sgc_catlases/{acc}_k31/cdbg.gxt",
        "outputs/metapangenome_sgc_catlases/{acc}_k31/bcalm.unitigs.db",
        "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/catlas.csv"
    resources: 
        mem_mb = 300000,
        time_min = 1440
    threads: 1
    conda: "envs/spacegraphcats.yml"
    params: outdir = "outputs/metapangenome_sgc_catlases"
    shell:'''
    python -m spacegraphcats build {input.conf} --outdir={params.outdir} --rerun-incomplete --nolock
    '''

#################################################################################
## Estimate abundance and perform dominating set differential abundance analysis
#################################################################################

rule metapangenome_spacegraphcats_catlas_cdbg_to_pieces_map:
    input:
        cdbg = "outputs/metapangenome_sgc_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/catlas.csv"
    output: "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/cdbg_to_pieces.csv"
    conda: "envs/spacegraphcats.yml"
    resources: 
        mem_mb = 16000,
        time_min = 720
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31_r10", 
    shell:'''
    scripts/cdbg_to_pieces.py {params.cdbg_dir} {params.catlas_dir}
    '''

rule tmp_cp_sgc_nbhds_w_lib_prefix:
    """
    the way that the abundance estimation code is written, it automatically
    1. performs across multiple samples to save on time with loading the catlas over and over
    2. automatically names files from the prefix of the basename of the input reads being estimated.
    The genome query neighborhoods all have the same prefix (the genome accession number) which would cause all of the files to overwrite one another.
    This rule temporarily copies those files to a new name, making the sample name the prefix.
    I could go back to the original metagenome and estimate abundances from those reads, but that takes a lot more time than just using query neighborhoods.
    I've tried to think of ways to fix this in the spacegraphcats code, but there wasn't an immediate answer for a better way to do this. 
    https://github.com/spacegraphcats/spacegraphcats/issues/451#issuecomment-987043304
    """
    input: "outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    output: temp("outputs/mgx_sgc_genome_queries_tmp/{acc}/{sample}.reads.gz")
    resources: 
        mem_mb = 500,
        time_min = 30
    threads: 1
    shell:'''
    cp {input} {output}
    '''
    
rule spacegraphcats_pangenome_catlas_estimate_abundances:
    input:
        cdbg = "outputs/metapangenome_sgc_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/catlas.csv",
        reads = expand("outputs/mgx_sgc_genome_queries_tmp/{{acc}}/{sample}.reads.gz", sample = SAMPLES)
    output: expand("outputs/metapangenome_sgc_catlases/{{acc}}_k31_r10_abund/{sample}.reads.gz.dom_abund.csv", sample = SAMPLES)
    conda: "envs/spacegraphcats.yml"
    resources: 
        mem_mb = 10000,
        time_min = 440
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31_r10", 
        out_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31_r10_abund", 
    shell:'''
    python -m spacegraphcats.search.count_dominator_abundance {params.cdbg_dir} {params.catlas_dir} --outdir {params.out_dir} {input.reads}
    '''

rule format_spacegraphcats_pangenome_catlas_abundances:
    input: 
        dom_abund = expand("outputs/metapangenome_sgc_catlases/{{acc}}_k31_r10_abund/{sample}.reads.gz.dom_abund.csv", sample = SAMPLES)
    output: 
        dom_abund="outputs/metapangenome_sgc_catlases/{acc}_k31_r10_abund/all_dom_abund.tsv",
        dom_info="outputs/metapangenome_sgc_catlases/{acc}_k31_r10_abund/dom_info.tsv",
        dom_abund_pruned="outputs/metapangenome_sgc_catlases/{acc}_k31_r10_abund/all_dom_abund_pruned.tsv"
    conda: "envs/tidyverse.yml"
    resources: 
        mem_mb = 200000,
        time_min = 440
    threads: 1
    script: "scripts/format_metapangenome_catlas_dom_abund.R"

rule corncob_for_dominating_set_differential_abund:
    input: 
        dom_abund_pruned="outputs/metapangenome_sgc_catlases/{acc}_k31_r10_abund/all_dom_abund_pruned.tsv",
        ntcard="outputs/mgx_ntcard/all_kmer_count.tsv",
        info = "inputs/metadata.csv"
    output: 
        all_ccs = "outputs/metapangenome_sgc_catlases_corncob/{acc}_all_ccs.tsv",
        sig_ccs = "outputs/metapangenome_sgc_catlases_corncob/{acc}_sig_ccs.tsv"
    resources: 
        mem_mb = 16000,
        time_min = 1440
    threads: 1
    conda: 'envs/corncob.yml'
    script: "scripts/corncob_dda.R"
