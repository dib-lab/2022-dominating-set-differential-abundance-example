import csv
import re
import pandas as pd

query_genomes = pd.read_csv("outputs/query_genomes_from_sourmash_gather/query_genomes.csv", header = 0)
GTDB_SPECIES = query_genomes['species'].tolist()
GTDB_SPECIES = [re.sub(' ', '_', i) for i in GTDB_SPECIES]
ACC = query_genomes['accession'].tolist()

# make a variable that binds GTDB species to the accession that the queries were seeded with;
# catlases are named after the accession, so these need to be bound so that each catlas is only annotated by the correct species
# Bound by two dashes ("--") so snakemake can parse wildcards correctly.
zip_to_list_of_lists = lambda ACC, GTDB_SPECIES: [list(tmp) for tmp in zip(ACC, GTDB_SPECIES)]
acc_to_gtdb_tmp = zip_to_list_of_lists(ACC, GTDB_SPECIES)

acc_to_gtdb = []
for i in acc_to_gtdb_tmp:
    acc_to_gtdb_joined = "--".join(i)
    acc_to_gtdb.append(acc_to_gtdb_joined)

ACC_SPECIES = acc_to_gtdb

class Checkpoint_GrabAnnotationAccessions:
    """
    Define a class to simplify file specification from checkpoint
    (e.g. solve for {acc} wildcard without needing to specify a function
    for each arm of the DAG that uses the acc wildcard).
    This approach is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self, gtdb_species):
        acc_csv = f'outputs/query_genomes_to_annotation_genomes/{gtdb_species}_gtdb_genomes.csv'
        assert os.path.exists(acc_csv)

        genome_accs = []
        with open(acc_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['name']
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {acc_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of rule 'grab_annotation_genome_accessions';
        # this will trigger exception until that rule has been run.
        checkpoints.grab_annotation_genome_accessions.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(w.gtdb_species)

        p = expand(self.pattern, acc=genome_accs, **w)
        return p

rule all:
    input: expand("outputs/metapangenome_sgc_catlases_corncob_annotations/{acc_species}_sig_ccs_annotated_no_cdbg.tsv", acc_species = ACC_SPECIES)

checkpoint grab_annotation_genome_accessions:
    input:
        lineages="inputs/gtdb-rs207.taxonomy.csv",
        query_genomes="outputs/query_genomes_from_sourmash_gather/query_genomes.csv"
    output: csv="outputs/query_genomes_to_annotation_genomes/{gtdb_species}_gtdb_genomes.csv",
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb = 4000,
        time_min = 30
    threads: 1
    script: "scripts/grab_annotation_genome_accessions.R"

rule make_annotation_genome_info_csv:
    output:
        csvfile = 'outputs/annotation_genomes/{acc}.info.csv'
    conda: "envs/genbank_genomes.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell: """
        python scripts/genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
    """
    
rule download_annotation_genome:
    input:
        csvfile = ancient('outputs/annotation_genomes/{acc}.info.csv')
    output:
        genome = "outputs/annotation_genomes/{acc}_genomic.fna.gz"
    resources:
        mem_mb = 500
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

rule bakta_download_db:
    output: "inputs/bakta_db/db/version.json"
    threads: 1
    resources:
        mem_mb = 4000,
        time_min = 240
    params: outdir = "inputs/bakta_db"
    conda: "envs/bakta.yml"
    shell:'''
    bakta_db download --output {params.outdir}
    '''

rule bakta_annotate_gtdb_genomes:
    input:
        fna=ancient("outputs/annotation_genomes/{acc}_genomic.fna.gz"),
        db="inputs/bakta_db/db/version.json",
    output:
        "outputs/annotation_genomes_bakta/{gtdb_species}/{acc}.faa",
        "outputs/annotation_genomes_bakta/{gtdb_species}/{acc}.gff3",
        "outputs/annotation_genomes_bakta/{gtdb_species}/{acc}.fna",
        "outputs/annotation_genomes_bakta/{gtdb_species}/{acc}.ffn",
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        time_min = 40
    conda: 'envs/bakta.yml'
    params:
        dbdir="inputs/bakta_db/db/",
        outdir_b = lambda wildcards: 'outputs/annotation_genomes_bakta/' + wildcards.gtdb_species,
    threads: 1
    shell:'''
    bakta --threads {threads} --db {params.dbdir} --prefix {wildcards.acc} --output {params.outdir_b} \
        --locus {wildcards.acc} --locus-tag {wildcards.acc} --keep-contig-headers {input.fna}
    '''

rule cat_annotated_sequences:
    input: ancient(Checkpoint_GrabAnnotationAccessions("outputs/annotation_genomes_bakta/{{gtdb_species}}/{acc}.ffn"))
    output: "outputs/annotation_genomes_annotated_combined/{gtdb_species}_all_annotated_seqs.fa"
    threads: 1
    resources:
        mem_mb=2000,
        time_min = 20
    shell:'''
    cat {input} > {output}
    '''

rule cluster_annotated_sequences:
    input: "outputs/annotation_genomes_annotated_combined/{gtdb_species}_all_annotated_seqs.fa"
    output: "outputs/annotation_genomes_annotated_clustered/{gtdb_species}_clustered_annotated_seqs.fa"
    threads: 1
    resources:
        mem_mb=16000,
        time_min=30
    conda: "envs/cdhit.yml"
    shell:'''
    cd-hit-est -c .95 -d 0 -i {input} -o {output}
    '''
    
rule translate_clustered_sequences_for_annotations:
    input: "outputs/annotation_genomes_annotated_clustered/{gtdb_species}_clustered_annotated_seqs.fa"
    output: "outputs/annotation_genomes_annotated_clustered/{gtdb_species}_clustered_annotated_seqs.faa"
    conda: 'envs/emboss.yml'
    resources:
        mem_mb = 16000,
        time_min=20
    threads: 2
    shell:'''
    transeq {input} {output}
    '''

rule eggnog_download_db:
    output: "inputs/eggnog_db/eggnog.db"
    threads: 1
    resources:
        mem_mb = 4000,
        time_min=2880
    params: datadir = "inputs/eggnog_db"
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir {params.datadir}
    '''

rule eggnog_annotate_clustered_sequences:
    input:
        faa = "outputs/annotation_genomes_annotated_clustered/{gtdb_species}_clustered_annotated_seqs.faa",
        db = 'inputs/eggnog_db/eggnog.db'
    output: "outputs/annotation_genomes_annotated_clustered_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations"
    conda: 'envs/eggnog.yml'
    resources:
        mem_mb = 32000,
        time_min=2880
    threads: 8
    params:
        outdir_e = lambda wildcards: "outputs/annotation_genomes_annotated_clustered_eggnog/" + wildcards.gtdb_species,
        dbdir = "inputs/eggnog_db"
    shell:'''
    mkdir -p tmp/
    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.gtdb_species} \
       --output_dir {params.outdir_e} -m hmmer -d none --tax_scope auto \
       --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       -d 2 --data_dir {params.dbdir}
    '''

rule sketch_metapangenome_reference:
    input: "outputs/annotation_genomes_annotated_clustered/{gtdb_species}_clustered_annotated_seqs.fa",
    output: "outputs/annotation_genomes_annotated_clustered_sigs/{gtdb_species}_clustered_annotated_seqs.sig"
    resources:
        mem_mb = 500,
        time_min = 60
    threads: 1
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=31,scaled=1000,abund -o {output} {input}
    '''

rule make_metapangenome_sgc_multifasta_conf_files:
    input:
        reads = "outputs/mgx_sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
        ref_genes = expand("outputs/annotation_genomes_annotated_clustered/{gtdb_species}_clustered_annotated_seqs.fa", gtdb_species = GTDB_SPECIES),
        ref_sig = expand("outputs/annotation_genomes_annotated_clustered_sigs/{gtdb_species}_clustered_annotated_seqs.sig", gtdb_species = GTDB_SPECIES)
    output:
        conf = "outputs/sgc_conf/{acc}_r10_multifasta_conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        row = query_genomes.loc[query_genomes['accession'] == wildcards.acc]
        gtdb_species_tmp = row['species'].values[0]
        gtdb_species = re.sub(" ", "_", gtdb_species_tmp)
        for gtdb_species_wc in wildcards.gtdb_species:
            if gtdb_species == gtdb_species_wc:
                reads = "outputs/mgx_sgc_genome_queries_hardtrim/" + wildcards.acc + ".hardtrim.fa.gz"
                ref_genes = "outputs/annotation_genomes_annotated_clustered/" + gtdb_species + "_clustered_annotated_seqs.fa"
                ref_sig =  "outputs/annotation_genomes_annotated_clustered_sigs/" + gtdb_species + "_clustered_annotated_seqs.sig"
                with open(output.conf, 'wt') as fp:
                    print(f"""\
catlas_base: {wildcards.acc}
input_sequences:
- {reads}
radius: 10
paired_reads: true
multifasta_reference:
- {ref_genes}
multifasta_scaled: 1000
multifasta_query_sig: {ref_sig}
""", file=fp)


## ADD A RULE THAT TOUCHES THE EXISTING CATLAS FIRST??
rule spacegraphcats_metapangenome_catlas_multifasta_annotate:
    input:
        conf = "outputs/sgc_conf/{acc}_r10_multifasta_conf.yml",
        catlas = "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/catlas.csv"
    output:
        annot="outputs/metapangenome_sgc_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_annot.csv",
        record="outputs/metapangenome_sgc_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_by_record.csv",
    params:
        outdir = "outputs/metapangenome_sgc_catlases/",
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 32000,
        time_min = 720
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} multifasta_query --nolock --outdir {params.outdir} --rerun-incomplete
    '''

rule spacegraphcats_metapangenome_catlas_cdbg_to_pieces_map:
    input:
        cdbg = "outputs/metapangenome_sgc_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/catlas.csv"
    output: "outputs/metapangenome_sgc_catlases/{acc}_k31_r10/cdbg_to_pieces.csv"
    conda: "envs/spacegraphcats.yml"
    resources: 
        mem_mb = 16000,
        time_min = 440
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/metapangenome_sgc_catlases/" + wildcards.acc + "_k31_r10", 
    shell:'''
    scripts/cdbg_to_pieces.py {params.cdbg_dir} {params.catlas_dir}
    '''

rule join_annotations_to_dominating_set_differential_abundance_analysis_results:
    input:
        cdbg_to_pieces="outputs/metapangenome_sgc_catlases/{acc}_k31_r10/cdbg_to_pieces.csv",
        eggnog = "outputs/annotation_genomes_annotated_clustered_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations", 
        cdbg_annot = "outputs/metapangenome_sgc_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_annot.csv",
        sig_ccs = "outputs/metapangenome_sgc_catlases_corncob/{acc}_sig_ccs.tsv",
        dom_info = "outputs/metapangenome_sgc_catlases/{acc}_k31_r10_abund/dom_info.tsv"
    output:
        sig_ccs_annot = "outputs/metapangenome_sgc_catlases_corncob_annotations/{acc}--{gtdb_species}_sig_ccs_annotated_no_cdbg.tsv"
    threads: 1
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb=8000,
        time_min = 60
    script: "scripts/join_annotations_to_dominating_set_differential_abundance_analysis_results.R"
