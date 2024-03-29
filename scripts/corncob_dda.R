library(magrittr)
library(dplyr)
library(tibble)
library(readr)
library(corncob)

# run corncob on dominating set pieces to determine which portions of 
# the pangenome graph are differentially abundant in IBD subtypes. 

# From old script, may no longer be relevant:
# Note that this script back-adds
# columns for samples that had no reads map (and therefore have empty files)
# during salmon quantification. 

# read in metadata --------------------------------------------------------

# read in sample metadata
info <- read_csv(snakemake@input[["info"]]) %>%
  select(sample, var) %>%
  distinct() 

# read in abundtrim lib size info
libsizes <- read_tsv(snakemake@input[['ntcard']])

# join metadata
info <- left_join(info, libsizes, by = "sample")

# import counts -----------------------------------------------------------

count_info <- read_tsv(snakemake@input[['dom_abund_pruned']])
count_info <- as.data.frame(count_info)
rownames(count_info) <- count_info$dom_id
count_info <- as.data.frame(count_info[ , -1])

# format counts for bdml ----------------------------------------------------

# transpose counts
count_info_t <- count_info %>% 
  as.data.frame %>% 
  t %>% 
  as_tibble

# re-add sample names
count_info_t %<>%
  add_column("sample" = colnames(count_info), .before=TRUE)

## Join up the data -- observations are rows and variables are columns
# join takes too long, use a different method
info <- info[order(match(info$sample, count_info_t$sample)), ]
# check that order matches; fail if not
stopifnot(all.equal(info$sample, count_info_t$sample))
df <- as.data.frame(cbind(info, count_info_t[ , -1])) 
# change levels of diagnosis so nonIBD is default
df$var <- factor(df$var, levels = c("nonIBD", "CD")) # THIS NEEDS TO BE UPDATE TO WHATEVER THE USER IS DOING

# Run corncob -------------------------------------------------------------

## function to test for differences in all AA seqs
## note inherits df from global env
fit_corncob <- function(col_num) {
  # record protein being queried
  domset_piece_name <- df %>%
    select(all_of(col_num)) %>%
    colnames()
  # capture the warning message output by corncob when a group contains zero counts
  tc <- textConnection("messages","w")
  sink(tc, type="message")
  # run corncob with LRT
  corncob_out <- df %>%
    select(num_kmers, var, all_of(col_num)) %>%
    rename(ww = 3) %>%
    corncob::bbdml(formula = cbind(ww, num_kmers - ww) ~ var,
                   #formula_null = cbind(ww, num_kmers - ww) ~ variable, # use this to block if you have additional variables
                   phi.formula = ~ 1,
                   phi.formula_null = ~ 1,
                   data = .,
                   test = "LRT", 
                   boot = FALSE) %>%
    summary()
  # close message text connection
  sink(NULL, type="message")
  close(tc)
  # make a dataframe for the final results
  corncob_coeff <- data.frame(mu = rownames(corncob_out$coefficients), 
                              estimate = corncob_out$coefficients[ , 'Estimate'],
                              standard_error = corncob_out$coefficients[ , 'Std. Error'],
                              t_value = corncob_out$coefficients[ , "t value"],
                              p_value = corncob_out$coefficients[ , "Pr(>|t|)"],
                              domset_piece = domset_piece_name,
                              separation_in_abund_model = ifelse(length(messages) == 0, "none", "separation"))
  return(corncob_coeff)
}

# run function on all domset pieces
all_ccs <- sapply(4:ncol(df), fit_corncob, simplify=F)
all_ccs <- do.call(rbind, all_ccs)
write_tsv(all_ccs, path = snakemake@output[["all_ccs"]])

# perfom pvalue adustment
sig <- all_ccs %>%
  filter(!grepl("Intercept", mu)) %>% # remove intercept before performing pvalue correction
  group_by(mu) %>%
  mutate(bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
  filter(bonferroni < .05)
write_tsv(sig, path = snakemake@output[["sig_ccs"]])
