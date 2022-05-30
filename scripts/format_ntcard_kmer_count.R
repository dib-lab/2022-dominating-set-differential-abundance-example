library(dplyr)
library(readr)
library(purrr)

files <- unlist(snakemake@input[["fstat"]])
fstat <- files %>%
  set_names() %>%
  map_dfr(read_tsv, n_max = 1, col_names = c("k", "F", "num_kmers"), 
          .id = "sample") %>%
  mutate(sample = gsub("\\.fstat", "", basename(sample))) %>%
  select(sample, num_kmers)

write_tsv(fstat, snakemake@output[["tsv"]])
