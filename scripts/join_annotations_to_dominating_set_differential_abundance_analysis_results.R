library(dplyr)
library(readr)
library(tidyr)

eggnog <- read_tsv(snakemake@input[['eggnog']], skip = 4, show_col_types = F) %>%
  rename(query_name = `#query`) %>%
  mutate(query_name = gsub("_1$", "", query_name)) # remove trailling one added by eggnog or clustering

cdbg_annot <- read_csv(snakemake@input[['cdbg_annot']], show_col_types= F) %>%
  separate(record_name, into = c("record_id", "tmp"), remove = F, sep = " ") %>%
  select(-tmp)

sig_ccs <- read_tsv(snakemake@input[['sig_ccs']]) %>%
  rename(dom_id = aa_seq)

cdbg_to_pieces <- read_csv(snakemake@input[['cdbg_to_pieces']])
dom_info <- read_tsv(snakemake@input[['dom_info']])

# filter to sig_ccs
sig_ccs <- left_join(sig_ccs, cdbg_to_pieces, by = c("dom_id" = "dominator"))
sig_ccs <- left_join(sig_ccs, cdbg_annot, by = c("cdbg_node" = "cdbg_id"))
sig_ccs <- left_join(sig_ccs, eggnog, by = c("record_id" = "query_name"))

# rm cdbg_node
sig_ccs <- sig_ccs %>%
  select(-cdbg_node) %>%
  distinct()

sig_ccs <- left_join(sig_ccs, dom_info, by = "dom_id")

write_tsv(sig_ccs, snakemake@output[['sig_ccs_annot']])
