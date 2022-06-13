library(readr)
library(dplyr)

query_genomes <- read_tsv(snakemake@input[['query_genomes']])
lineages <- read_csv(snakemake@input[['lineages']])

# grab the species of interest
accession_w <- unlist(snakemake@wildcards[['acc']])
query_genome_w <- query_genomes %>%
  filter(accessions == accession_w)

# write a tsv file that contains the lineage information for each set of species
# specified in the metadata file

lineages_w <- lineages %>%
  filter(species %in% query_genome_w$species) %>%
  select(name = ident, superkingdom, phylum, class, order, family, genus, species)

write_csv(lineages_w, snakemake@output[["csv"]])
