source("code/log_smk.R")
library(mikropml)

doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])
options(future.globals.maxSize = Inf)

data_raw <- readr::read_csv(snakemake@input[["csv"]])
dim(data_raw)

data_processed <- preprocess_data(data_raw,
                                  outcome_colname = snakemake@params[['outcome_colname']],
                                  group_neg_corr = TRUE,
                                  remove_var = 'nzv')
summary(data_processed)

saveRDS(data_processed, file = snakemake@output[["rds"]])
