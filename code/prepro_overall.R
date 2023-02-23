# Builds individual phenotype file and generates necessary directories for snakemake workflow
# This file would need to be updated based on the measures that the investigator is taking to adjust the continuous measure
# In this example, cytokines of interest are adjusted based on association with weighted elix score
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#FILES ----
pheno <- readr::read_csv(file = snakemake@input[["metadata"]])
phenotype_cols <- snakemake@wildcards[["phenotype"]]

paths <- paste0(c("data/pheno/",
                  "results/",
                  "benchmarks/",
                  "figures/",
                  "log/",
                  "data/mikropml/"),
                phenotype_cols)

for(i in 1:length(paths)){
  
  if(dir.exists(paths[i]) == FALSE){
    
    dir.create(paths[i])
    
  }else{
    
    print(paste0("directory ", paths[i], " already exists"))
    
  }
}

new_dir3 <- paste0("results/", phenotype_cols, "/runs")

if(dir.exists(new_dir3) == FALSE){
  
  dir.create(new_dir3)
    
  }else{
    
    print(paste0("directory ", new_dir3, " already exists"))
    
  }

pheno_sub <- pheno %>%
  select(genome_id,
         all_of(phenotype_cols)) %>%
  drop_na()
  
write_tsv(pheno_sub,
          file = paste0("data/pheno/", phenotype_cols, "/handsoff.tsv"))
