source("code/log_smk.R")

doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])
options(future.globals.maxSize = Inf)

data_processed <- readRDS(snakemake@input[["rds"]])$dat_transformed

default_params <- mikropml::get_hyperparams_list(dataset = data_processed,
                                                 method = snakemake@params[["method"]])

lambda_alt <- c(0.00001, 0.0001, 0.001, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.1, 0.5, 1, 10)

alpha_var <- snakemake@params[["alpha_var"]]
lambda_var<- snakemake@params[["lambda_var"]]

alpha_assign <- if(alpha_var == "default"){
  default_params$alpha
}else{
  seq(0, 1, as.numeric(alpha_var))
}

lambda_assign <- if(lambda_var == "default"){
  default_params$lambda
}else{
  lambda_alt
}

new_hp <- list("alpha" = alpha_assign,
               "lambda" = lambda_assign)

ml_results <- mikropml::run_ml(
  dataset = data_processed,
  method = snakemake@params[["method"]],
  outcome_colname = snakemake@params[['outcome_colname']],
  find_feature_importance = TRUE,
  hyperparameters = new_hp,
  kfold = as.numeric(snakemake@params[['kfold']]),
  seed = as.numeric(snakemake@params[["seed"]])
)

saveRDS(ml_results$trained_model, file = snakemake@output[["model"]])
readr::write_csv(ml_results$performance, snakemake@output[["perf"]])
readr::write_csv(ml_results$feature_importance, snakemake@output[["feat"]])
