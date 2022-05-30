library(dplyr)
library(readr)
library(purrr)
library(tidyr)

# read in gather results and gtdb lineages --------------------------------

gtdb_lineages <- read_csv(snakemake@input[['lineages']])
#gtdb_lineages <- read_csv("~/github/2021-metapangenome-example/inputs/gtdb-rs202.taxonomy.v2.csv")

gather_results <- unlist(snakemake@input[['gather']]) %>%
#gather_results <- Sys.glob("~/github/2021-metapangenome-example/outputs/sample_gather/*genomic.csv") %>%
  set_names() %>%
  map_dfr(read_csv, col_types = c("dddddlllcccddddcccd"), .id = "sample") %>%
  mutate(sample = gsub("_gather_gtdb-rs202-genomic.csv", "", basename(sample))) %>%
  separate(col = name, into = c("accession"), remove = F, sep = " ", extra = "drop")

gather_results <- left_join(gather_results, gtdb_lineages, by = c("accession" = "ident"))


# select query genomes ----------------------------------------------------

min_sample_frac <- snakemake@params[['min_sample_frac']]

# count the number of samples each genome was detected in
gather_results_genome_tally <- gather_results %>%
  group_by(name) %>%
  tally() %>%
  select(name, n_samples_genome_detected = n)

# filter to genomes that were detected in min_sample_frac samples
query_genome_results <- gather_results %>%
  left_join(gather_results_genome_tally, by = "name") %>%
  mutate(n_samples = length(unique(sample)),
         sample_frac =  n_samples_genome_detected / n_samples) %>%
  filter(sample_frac >= min_sample_frac) %>%
  select(accession, superkingdom, phylum, class, order, family, genus, species) %>%
  distinct()

write_csv(query_genomes_results, snakemake@output[['query_genomes']])
