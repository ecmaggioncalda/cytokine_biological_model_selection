#mikropml input file generation
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#PHENO ----
paste0(snakemake@input[["pheno"]])
pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t",
                           col_names = TRUE)

pheno_merge <- as.data.frame(pheno)
rownames(pheno_merge) <- pheno_merge$genome_id
pheno_merge <- pheno_merge[, -1, drop = FALSE]

print("pheno_merge has been successfully created")
dim(pheno_merge)

#GENO ----
paste0(snakemake@params[['path']])
geno <- readr::read_delim(file = snakemake@params[['path']],
                          col_names = TRUE)

print("geno matrix has been successfully read in")
dim(geno)

#FEATURE SELECTION ----
feature <- read_csv(file = snakemake@input[["feature_file"]])

feature_sub <- feature %>%
  select(full_names,
         all_of(colnames(pheno_merge))) %>%
  filter(if_any(where(is.numeric), ~.x == 1))

print("feature matrix has been successfully read in and subsetted for phenotype")
dim(feature_sub)

index <- sapply(feature_sub$full_names,
                function(x){
                  
                  which(x == geno$variant)
                  
                })

if(length(index) != length(feature_sub$full_names)){
  stop("not all indexes were identified in the geno matrix")
}
if(any(is.na(index))){
  stop("not all indices were identified in the geno matrix")
}

geno_sub <- geno[index,] %>%
  select(variant,
         all_of(rownames(pheno_merge))) %>%
  mutate(variant = gsub("_$", "", variant)) %>%
  column_to_rownames("variant")

geno_merge <- t(geno_sub)

if(sum(rownames(pheno_merge) %in% rownames(geno_merge)) != length(rownames(pheno_merge))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(pheno_merge), rownames(geno_merge))

geno_ordered <- geno_merge[index, , drop = FALSE]

if(sum(rownames(pheno_merge) == rownames(geno_ordered)) != length(rownames(pheno_merge))){
  stop("mismatch between pheno and geno contents")
}

complete_frame <- cbind(pheno_merge,
                        geno_ordered)

print("complete frame generated with pheno:geno, export to csv for mikropml preprocessing")

#GENERATE FILES ----
#patient and genome factors
write_csv(complete_frame,
          file = snakemake@output[['file_name']],
          col_names = TRUE)
