import pandas as pd
import glob
import os

metadata_file = config["metadata_file"]
m = pd.read_csv(metadata_file, header = 0)
SAMPLES = m['sample'].unique().tolist()


########################################
## PREPROCESSING
########################################

rule mgx_fastp_reads:
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
   This rule assumes a human host for a metagenome.
   If your metagenome did not come from a human, it's probably still a good idea to try and subtract the host genome.
   While here we used a masked version of the human genome, using the genome itself as the host will still work (no promises around mitochondira/chlorophyll).
   Replace the path to the host file with the path to your host genome file.
   We've provided download scripts for a variety of common mammalian hosts in host_sequences.snakefile (mostly because we had then laying around); any genome sequence could be used here however.
   These sequences can be used as drop-in replacements for the "host" sequence indicated below. 
   """
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
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
    input: 
        'outputs/mgx_bbduk/{sample}_R1.nohost.fq.gz',
        'outputs/mgx_bbduk/{sample}_R2.nohost.fq.gz'
    output: "outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources: 
        mem_mb = 64000,
        time_min = 720
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

##########################################################
## Determine taxonomic profile of metagenomes &
## identified species that are present in all metagenomes
##########################################################

rule mgx_sourmash_sketch:
    input: "outputs/mgx_abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/mgx_sigs/{sample}.sig"
    params: "outputs/mgx_sigs/{sample}.sig"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        time_min=1200,
    conda: "envs/sourmash.yml"
    shell:"""
        sourmash sketch dna -p k=21,31,51 scaled=2000,abund -o {output} --name {wildcards.sample} {input}
    """

rule download_sourmash_gather_database:
    output: "inputs/gtdb-rs207.genomic-reps.k31.zip"
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
    input:
        sig="outputs/sample_sigs/{sample}.sig",
        db="inputs/gtdb-rs207.genomic-reps.dna.k31.zip",
    output: 
        csv="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv",
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 128000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_gather.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

checkpoint mgx_select_species_shared_across_samples:
    """
    
    """
    input:
    output:
        lineages="",

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
    input: csvfile = ancient('outputs/genbank_genomes/{acc}.info.csv')
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

checkpoint query_genomes_charcoal_decontaminate:
    """
    Because a relatively low containment index (0.01) between a query and a metagenome is necessary to recover ~20-40% of an organisms genome from the metagenome, contamination could potentially lead to recovery of many off target sequences. 
    Charcoal decontaminates the query genomes prior to running the assembly graph queries.
    """
    input:
        genomes = ancient(Checkpoint_GatherResults("outputs/query_genomes/{acc}_genomic.fna.gz")),
        genome_list = "outputs/charcoal_conf/charcoal.genome-list.txt",
        conf = "inputs/charcoal-conf.yml",
        genome_lineages="inputs/gtdb-rs207.taxonomy.csv.gz"
        db="inputs/gtdb-rs207.genomic-reps.k31.zip"
        db_lineages="inputs/gtdb-rs207.taxonomy.csv.gz"
    output: directory("outputs/query_genomes_charcoal/")
        #"outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz"
    resources:
        mem_mb = 64000
    threads: 8
    conda: "envs/charcoal.yml"
    shell:'''
    python -m charcoal run {input.conf} -j {threads} clean --nolock --latency-wait 15 --rerun-incomplete
    '''

def checkpoint_query_genomes_charcoal_decontaminate:
    # Expand checkpoint to get query genome accs, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.query_genomes_charcoal_decontaminate.get(**wildcards).output[0]
    file_names = expand("outputs/query_genomes_charcoal/{acc}_genomic.fna.gz.clean.fa.gz",
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz")).acc)
    return file_names

rule mgx_make_sgc_conf:

rule mgx_spacegraphcats_build_catlas:
    """
    The RAM necessary to run this rule is mostly tied to the construction of the cDBG. 
    This becomes more RAM-intensive when there are more k-mers in a sample.
    As such, the RAM necessary to run this step increases with the complexity/taxonomic diversity of the sample.
    We include some suggestions based on using this approach with metagenomes from a variety of environments.
    RAM suggestions:
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
    output:
        reads = directory("outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/")
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

##########################################################
## Build metapangenome graphs
##########################################################

rule diginorm_spacegraphcats_query_genomes:
    input: expand("outputs/mgx_sgc_genome_queries/{sample}_k31_r1_search_oh0/{{acc}}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz", sample = SAMPLES)
    output: "outputs/mgx_sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    resources:
        mem_mb = 164000,
        time_min = 1440
    threads: 1
    conda: "envs/env.yml"
    shell:'''
    zcat {input} | normalize-by-median.py -k 20 -C 20 -M 164e9 --gzip -o {output} -
    '''

rule hardtrim_spacegraphcats_query_genomes:
    input: "outputs/mgx_sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    output: "outputs/mgx_sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz"
    resources:
        mem_mb = 24000,
        time_min = 1440
    threads: 1
    conda: "envs/env.yml"
    shell:'''
    trim-low-abund.py -C 4 -M 20e9 -k 31 {input} --gzip -o {output}
    '''

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
    input:
        reads = "outputs/mgx_sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
        conf = "outputs/sgc_conf/{acc}_r10_conf.yml"
    output: 
        "outputs/metapangenome_sgc_catlases/{acc}_k31/cdbg.gxt",
        "outputs/metapangenome_sgc_catlases/{acc}_k31/bcalm.unitigs.db",
        "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/catlas.csv"
    resources: mem_mb = 300000
    conda: "envs/spacegraphcats.yml"
    params: outdir = "outputs/sgc_pangenome_catlases"
    shell:'''
    python -m spacegraphcats build {input.conf} --outdir={params.outdir} --rerun-incomplete --nolock
    '''

#################################################################################
## Estimate abundance and perform dominating set differential abundance analysis
#################################################################################

rule spacegraphcats_pangenome_catlas_cdbg_to_pieces_map:
    input:
        cdbg = "outputs/sgc_pangenome_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
    output: "outputs/sgc_pangenome_catlases/{acc}_k31_r10/cdbg_to_pieces.csv"
    conda: "envs/spacegraphcats2.yml"
    resources: 
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10", 
    shell:'''
    scripts/cdbg_to_pieces.py {params.cdbg_dir} {params.catlas_dir}
    '''

rule tmp_cp_sgc_nbhds_w_lib_prefix:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    output: temp("outputs/sgc_genome_queries_tmp/{acc}/{library}.reads.gz")
    resources: 
        mem_mb = 500,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    cp {input} {output}
    '''
    
# TR TODO: update env to PR 303, or update sgc latest if merged. Since dom_abund is checked out, this might work like this...
rule spacegraphcats_pangenome_catlas_estimate_abundances:
    input:
        cdbg = "outputs/sgc_pangenome_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv",
        reads = expand("outputs/sgc_genome_queries_tmp/{{acc}}/{library}.reads.gz", library = LIBRARIES)
        #reads = expand("outputs/abundtrim/{library}.abundtrim.fq.gz", library = LIBRARIES)
    output: expand("outputs/sgc_pangenome_catlases/{{acc}}_k31_r10_abund/{library}.reads.gz.dom_abund.csv", library = LIBRARIES)
    conda: "envs/spacegraphcats_dom.yml"
    resources: 
        mem_mb = 10000,
        tmpdir = TMPDIR
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10", 
        out_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10_abund", 
    shell:'''
    /home/tereiter/github/spacegraphcats/scripts/count-dominator-abundance.py {params.cdbg_dir} {params.catlas_dir} --outdir {params.out_dir} {input.reads}
    '''

rule format_spacegraphcats_pangenome_catlas_abundances:
    input: 
        dom_abund = expand("outputs/sgc_pangenome_catlases/{{acc}}_k31_r10_abund/{library}.reads.gz.dom_abund.csv", library = LIBRARIES)
    output: 
        dom_abund="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/all_dom_abund.tsv",
        dom_info="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/dom_info.tsv",
        dom_abund_pruned="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/all_dom_abund_pruned.tsv"
    conda: "envs/tidy.yml"
    resources: 
        mem_mb = 200000,
        tmpdir = TMPDIR
    threads: 1
    script: "scripts/format_pangenome_catlas_dom_abund.R"

rule install_corncob:
    output: corncob = "outputs/sgc_pangenome_catlases_corncob/corncob_install.txt"
    resources: 
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    conda: 'envs/corncob.yml'
    script: "scripts/install_corncob.R"

rule corncob_for_dominating_set_differential_abund:
    input: 
        corncob="outputs/sgc_pangenome_catlases_corncob/corncob_install.txt",
        dom_abund_pruned="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/all_dom_abund_pruned.tsv",
        ntcard="outputs/ntcard/all_kmer_count.tsv",
        info = "inputs/working_metadata.tsv"
    output: 
        all_ccs = "outputs/sgc_pangenome_catlases_corncob/{acc}_all_ccs.tsv",
        sig_ccs = "outputs/sgc_pangenome_catlases_corncob/{acc}_sig_ccs.tsv"
    resources: 
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    conda: 'envs/corncob.yml'
    script: "scripts/corncob_dda.R"


#######################################################################
## Make reference multifasta sequence and annotate metapangenome graph
#######################################################################
