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
# only use intersection of samples between life hist and 16s data sets
life_hist_subset <- 
  life_hist %>%
  filter(Sample %in% asv_samples) %>%
  filter(!is.na(avg_final_size), !is.na(avg_cumulative_reproduction))

asv_wide <- asv_wide %>% filter(Sample %in% life_hist_subset$Sample)

# number of samples per treatment groups is ... weird bc of split-brood design
# could make problem with cross-validation if we use these as predictors
# maybe could use generations 1:5 from 15 gen expt as out of sample test?
life_hist_subset %>%
  group_by(generation, trt_type) %>% 
  summarise(n())

life_hist_subset %>%
  group_by(treatment) %>% 
  summarise(n())

# save each df arranged in alpha order by sample
# necessary because function for elastic net / lasso require vector reponse and matrix X (no data argument)
life_hist_subset <- life_hist_subset %>% arrange(Sample)
asv_wide <- asv_wide %>% arrange(Sample)

# check order of samples
all.equal(asv_wide$Sample, life_hist_subset$Sample)

# load package for elastic net regression
library(glmnet)

# make sparse design matrix
library(Matrix)
asv_mat <- as(as.matrix(asv_wide[,-1]), "sparseMatrix")

# 40 % of size of original matrix
# will speed up model fititing
# need to learn how to make sparse model matrix to accomodate interactions in the future
object.size(asv_mat) / object.size(asv_wide[,-1])

# set number of grid points for alpha
n_alpha <- 10

# init list to store cv objects
fit_list <- vector("list", length = length(0:10))

# set foldid
set.seed(123)
foldid <- sample(rep(1:10,length.out = nrow(life_hist_subset))) # control which obs go to which folds in CV

# fit for all values of alpha ranging from ridge to lasso
for (iter in 0:10) {
fit_list[[iter+1]] <- cv.glmnet(y = life_hist_subset$avg_final_size,
                                   x = asv_mat,
                                   type.measure = "mse",
                              foldid = foldid,
                                   alpha = (iter / 10))
print(iter)
}


# get the minimum mean c
all_cv <- data.frame(alpha = 0:10 / 10
           , min_cv_mean = sapply(seq_along(fit_list), function(fit) min(fit_list[[fit]]$cvm))
           #, cv_sd = sapply(seq_along(fit_list), function(fit) min(fit_list[[fit]]$cvsd))
           )

plot(min_cv_mean ~ alpha, data = all_cv)

best_mod <- fit_list[[2]]
best_mod

best_coef <- coef(best_mod, s = "lambda.1se")[-1] # drop intercept

best_coef = as.data.frame(as.matrix(best_coef))
indices = which(best_coef != 0)
key_otu = colnames(asv_mat)[indices]
length(key_otu)

coef_df <- data.frame(OTU = key_otu,
                      beta = best_coef[indices,])

asv <- asv %>% mutate(selected = OTU %in% key_otu)

asv_selected <- left_join(asv, coef_df)

# plot stuff
asv_selected %>%
  ggplot(aes(group = Order, y = beta)) +
  geom_boxplot()

# 
coef_df %>%
  arrange(beta) %>% hist

hist(coef_df$beta)

# find largest magnitude
asv_selected %>%
  filter(!is.na(beta)) %>%
  arrange(beta) %>%
  slice_head(n = 10) %>%
  select(Kingdom:beta)

asv_selected %>%
  filter(!is.na(beta)) %>%
  arrange(beta) %>%
  slice_tail(n = 10) %>%
  select(Kingdom:beta)
