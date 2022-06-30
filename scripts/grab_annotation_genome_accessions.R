library(readr)
library(dplyr)

query_genomes <- read_csv(snakemake@input[['query_genomes']])
lineages <- read_csv(snakemake@input[['lineages']])

# loop over species in query_genomes and write a csv file that contains
# all accessions for genomes of the same species in GTDB
for(species_w in query_genomes$species){
  # filter lineages to species of interest 
  lineages_w <- lineages %>%
    filter(species %in% species_w) %>%
    select(name = ident, superkingdom, phylum, class, order, family, genus, species)
  
  # edit species name to contain no space
  species_no_space_w <- gsub(" ", "_", species_w)
 
  # write to csv file
  write_csv(lineages_w, paste0("outputs/query_genomes_to_annotation_genomes/", species_no_space_w, "_gtdb_genomes.csv"))
}
