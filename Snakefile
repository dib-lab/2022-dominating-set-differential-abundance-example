import pandas as pd
import csv
import glob
import os

m = pd.read_csv("inputs/metadata.csv", header = 0)
SAMPLES = m['sample'].unique().tolist()

class Checkpoint_GrabAccessions:
    """
    Define a class to determine the query accessions to use based on sourmash gather results.
    This approach is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self):
        acc_csv = 'outputs/query_genomes_from_sourmash_gather/query_genomes.csv'
        assert os.path.exists(acc_csv)

        genome_accs = []
        with open(acc_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['accession']
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {acc_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of rule 'mgx_select_query_genomes_shared_across_samples'; 
        # this will trigger exception until that rule has been run.
        checkpoints.mgx_select_query_genomes_shared_across_samples.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(**w)

        p = expand(self.pattern, acc=genome_accs, **w)
        return p

rule all:
    input: Checkpoint_GrabAccessions("outputs/metapangenome_sgc_catlases_corncob/{acc}_sig_ccs.tsv")
        

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
        json = 'outputs/mgx_fastp/{sample}.json'
    conda: 'envs/fastp.yml'
    threads: 1
    resources: 
        mem_mb = 8000,
        time_min = 720
    shell:'''
    fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -q 4 -j {output.json} -l 31 -c
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

checkpoint mgx_select_query_genomes_shared_across_samples:
    """
    Read in the gather results for all samples, and select query genomes that are present in some minimum fraction of samples.
    The fraction of samples a genome must be detected in is set in params.min_sample_frac, and by default is set to 0.6 (it should probably be set to 1, but these samples are so dramatically reduced in size that 0.6 made more sense for the test data set).
    We intentionally use the reps database so that there will be more shared genomes detected between samples, at the expense of the more of the sample being identifiable.
    (E.g. if E. coli is present in all samples, we'll get the same E. coli reps match instead of the best matching genome, which may be different between samples.)
    """
    input:
        gather=expand("outputs/mgx_sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv", sample = SAMPLES),
        lineages="inputs/gtdb-rs207.taxonomy.csv.gz"
    output:
        query_genomes="outputs/query_genomes_from_sourmash_gather/query_genomes.csv",
    params: min_sample_frac = 0.6
    conda: 'envs/tidyverse.yml'
    resources:
        mem_mb = 4000,
        time_min = 60
    script: "scripts/mgx_select_query_genomes_shared_across_samples.R"
    
################################################################
##
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
        mem_mb=500
    shell:'''
    ls outputs/query_genomes/*gz | xargs -n 1 basename > {output} 
    '''

checkpoint query_genomes_charcoal_decontaminate:
    """
    Because a relatively low containment index (0.01) between a query and a metagenome is necessary to recover ~20-40% of an organisms genome from the metagenome, contamination could potentially lead to recovery of many off target sequences. 
    Charcoal decontaminates the query genomes prior to running the assembly graph queries.
    """
    input:
        genomes = ancient(Checkpoint_GrabAccessions("outputs/query_genomes/{acc}_genomic.fna.gz")), # expands the {acc} wildcard using the Checkpoint_GrabAccessions class
        genome_list = "outputs/charcoal_conf/charcoal.genome-list.txt",
        conf = "inputs/charcoal-conf.yml",
        genome_lineages="inputs/gtdb-rs207.taxonomy.csv.gz",
        db="inputs/gtdb-rs207.genomic-reps.dna.k31.zip",
        db_lineages="inputs/gtdb-rs207.taxonomy.csv.gz"
    output: directory("outputs/query_genomes_charcoal/")  # re-creates the {acc} wildcard, now assoc with charcoal output
        #"outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz"
    resources:
        mem_mb = 64000,
        time_min = 1440
    threads: 8
    conda: "envs/charcoal.yml"
    shell:'''
    python -m charcoal run {input.conf} -j {threads} clean --nolock --latency-wait 15 --rerun-incomplete
    '''

def checkpoint_query_genomes_charcoal_decontaminate(wildcards):
    # Expand checkpoint to get query genome accs, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.query_genomes_charcoal_decontaminate.get(**wildcards).output[0]
    file_names = expand("outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz",
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz")).acc)
    return file_names

rule mgx_make_sgc_conf:
    """
    Make configuration file that will be used for mgx_spacegraphcats_built_catlas and mgx_spacegraphcats_query_genomes_extract_reads
    """
    input: queries = checkpoint_query_genomes_charcoal_decontaminate # expands the {acc} wildcard from charcoal
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


rule mgx_spacegraphcats_build_catlas:
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
        reads = "outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz"
    output:
        "outputs/mgx_sgc_genome_queries/{sample}_k31_r10/catlas.csv",
        "outputs/mgx_sgc_genome_queries/{sample}_k31/cdbg.gxt",
        "outputs/mgx_sgc_genome_queries/{sample}_k31/bcalm.unitigs.db"
    params: outdir = "outputs/mgx_sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 32000,
        time_min = 1440
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} build --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

checkpoint mgx_spacegraphcats_query_genomes_extract_reads:
    input: 
        conf = ancient("outputs/sgc_conf/{sample}_k31_r1_conf.yml"),
        reads = "outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz"
    output: directory("outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/")  # re-creates the {acc} wildcard using sgc outputs
        #"outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz",
    params: outdir = "outputs/mgx_sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 10000,
        time_min = 1440
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

def checkpoint_mgx_spacegraphcats_query_genomes_extract_reads_1(wildcards):
    # Expand checkpoint to get query genome accs, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.mgx_spacegraphcats_query_genomes_extract_reads.get(**wildcards).output[0]
    #file_names = expand("outputs/mgx_sgc_genome_queries_hardtrim/{acc}.fa.gz",
    file_names = expand("outputs/mgx_sgc_genome_queries/{{sample}}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz",
                        #sample = wildcards.SAMPLES,
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz")).acc)
    return file_names

rule dummy_solve_sgc:
    input: checkpoint_mgx_spacegraphcats_query_genomes_extract_reads_1 # solve the SAMPLES/ACC wildcards from the checkpoint.
    #output: touch("outputs/mgx_spacegraphcats_query_genomes_extract_reads_{sample}_done.txt")
    output: touch("outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz") 

##########################################################
## Build metapangenome graphs
##########################################################

rule diginorm_spacegraphcats_query_genomes:
    input:
        #reads = checkpoint_mgx_spacegraphcats_query_genomes_extract_reads_1, # solve the SAMPLES/ACC wildcards from the checkpoint.
        #dummy = expand("outputs/mgx_spacegraphcats_query_genomes_extract_reads_{sample}_done.txt", sample = SAMPLES),
        reads = expand("outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/{{acc}}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz", sample = SAMPLES) # recreate ACC wildcard using Checkpoint_GrabAccessions to solve at end of workflow.
    output: "outputs/mgx_sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    resources:
        mem_mb = 164000,
        time_min = 1440
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    zcat {input.reads} | normalize-by-median.py -k 20 -C 20 -M 164e9 --gzip -o {output} -
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
    params: outdir = "outputs/sgc_pangenome_catlases"
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


#def checkpoint_mgx_spacegraphcats_query_genomes_extract_reads_2(wildcards):
#    # Expand checkpoint to get query genome accs, which will be used as queries for spacegraphcats extract
#    # checkpoint_output encodes the output dir from the checkpoint rule.
#    checkpoint_output = checkpoints.mgx_spacegraphcats_query_genomes_extract_reads.get(**wildcards).output[0]
#    file_names = expand("outputs/metapangenome_sgc_catlases_corncob/{acc}_sig_ccs.tsv",
#                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz")).acc)
#    return file_names


#######################################################################
## Make reference multifasta sequence and annotate metapangenome graph
#######################################################################
