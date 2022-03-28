library(phyloseq)
library(tidyverse)

# read in life history data and 16s data
asv <- readRDS("data/16s.rds") %>% transform_sample_counts(function(x) x/sum(x)) %>% psmelt()
life_hist <- read_csv("data/pooled_sample_life_history.csv")

# drop rows where no life hist is collected
# life_hist <- life_hist %>% filter(!is.na(avg_final_size), !is.na(avg_cumulative_reproduction))

# get names of samples in life history data frame
life_hist_samples <- life_hist %>% pull(Sample) 

# make asv dataframe wider: 
# 1 row for each sample
# 1 column for each OTU (not counting column w/ sample ID)
# elements of dataframe are relative abundance of each OTU in each sample
asv_wide <-
  asv %>% 
  select(Sample, OTU, Abundance) %>% 
  # filter(Sample %in% life_hist_samples) %>%
  pivot_wider(names_from = OTU, values_from = Abundance)

# 81 samples
# only 71 have life history info
dim(asv_wide)

# get labels of samples in asv
asv_samples <- asv_wide %>% pull(Sample)

# the asv data is missing these two samples that are present in life history data
missing_samples <- setdiff(life_hist_samples, asv_samples)
missing_samples

# which rows in life_hist are missing in 16s data?
life_hist[life_hist$Sample %in% missing_samples, ]

# this explains the discrepency between the number of rows in asv data and how many we would expect
# number of rows there should be
(n_OTU * nrow(life_hist))  # (# of unique OTU) * (# unique samples)

# number of rows in ASV data
nrow(asv) 

# it's because the asv data is missing two samples
# after accounting for these two missing samples we recover the number of rows in the asv dataframe
(n_OTU * (nrow(life_hist) - 2) ) == nrow(asv)

### Fit models
# only use subset of lifehistory that intersects with samples in 16s
life_hist_subset <- 
  life_hist %>%
  filter(Sample %in% asv_samples)

