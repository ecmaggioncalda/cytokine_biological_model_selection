#LIBRARIES ----
library(tidyverse)

#GWAS RESULTS ----
full_hits_prime <- read_csv("../../2023_01_16_combined_geno_summary/data/full_hits.csv")

treewas_tests <- unique(full_hits_prime$treewas.test)
hogwash_tests <- unique(full_hits_prime$hogwash.test)[!is.na(unique(full_hits_prime$hogwash.test))]
ml_tests <- unique(full_hits_prime$ml.test)
genos <- unique(full_hits_prime$geno)
cyt_groups <- unique(full_hits_prime$group)
cytokines <- unique(full_hits_prime$pheno)

full_hits_bool <- full_hits_prime %>%
  select(pheno,
         geno,
         group,
         names,
         hogwash.test,
         hogwash_Sig_Pval,
         treewas.test,
         treewas.sig) %>%
  distinct() %>%
  pivot_wider(names_from = hogwash.test,
              values_from = hogwash_Sig_Pval,
              values_fn = unique) %>%
  pivot_wider(names_from = treewas.test,
              values_from = treewas.sig,
              values_fn = unique) %>%
  mutate(`NA` = NULL) %>%
  distinct() %>%
  group_by(pheno,
           geno,
           group,
           names) %>%
  rowwise() %>%
  summarize(hogwash_sig = sum(c_across(c(all_of(hogwash_tests))), na.rm = TRUE),
            treewas_sig = sum(c_across(c(all_of(treewas_tests))), na.rm = TRUE)) %>%
  distinct() %>%
  mutate(sig_groupings = paste0(hogwash_sig, ".", treewas_sig))

full_hits <- full_hits_prime %>%
  full_join(full_hits_bool, by = intersect(colnames(full_hits_bool), colnames(full_hits_prime))) %>%
  distinct() %>%
  mutate(pheno.group = paste0(pheno, ".", group))

sig_hits <- full_hits %>%
  filter(sig_groupings != "0.0")

sig_hits_distinct <- sig_hits %>%
  select(pheno.group,
         geno,
         full_names,
         sig_groupings) %>%
  distinct() %>%
  mutate(full_names = gsub("_$", "", full_names))

sig_hits_full_names <- unique(sig_hits_distinct$full_names)

sig_hits_pivot <- sig_hits_distinct %>%
  select(pheno.group,
         full_names) %>%
  pivot_wider(names_from = pheno.group,
              values_from = pheno.group,
              values_fn = length,
              values_fill = 0)

write_csv(sig_hits_pivot,
          "data/minimal_filtered_features.csv")

#GENOS ----
core_df <- read_delim("../../data/core_mat_sift.tsv")
gene_df <- read_delim("../../data/gene_mat.tsv")
pan_df <- read_delim("../../data/pan_mat.tsv")
struct_df <- read_delim("../../data/pan_struct_mat.tsv")

geno_df <- bind_rows(core_df,
                     gene_df,
                     pan_df,
                     struct_df)

index <- sapply(sig_hits_full_names,
                function(x){
                  
                  which(x == geno_df$variant)
                  
                })

length(index) == length(sig_hits_full_names) #needs to be TRUE
any(is.na(index)) #needs to be FALSE

geno_df_sub <- geno_df[index,]

write_delim(geno_df_sub,
            "data/combined_mat.tsv")

#PHENOTYPE ----
pheno_dirs <- list.files("../../2023_01_04_snakemake_sift_core_analysis/cytokine_core_sift/data/pheno",
                         full.names = TRUE)
pheno_path <- unlist(lapply(pheno_dirs, function(x){list.files(x,
                                                               pattern = "*.tsv",
                                                               full.names = TRUE)}))
pheno <- read_delim(pheno_path[1])
colnames(pheno)[2] <- gsub("\\.tsv",
                           "",
                           paste0(str_split(pheno_path[1],
                                            "/",
                                            simplify = TRUE)[,7],
                                  ".",
                                  str_split(pheno_path[1],
                                            "/",
                                            simplify = TRUE)[,8]))

for(i in 2:length(pheno_path)){
  
  a <- read_delim(pheno_path[i])
  colnames(a)[2]<- gsub("\\.tsv",
                        "",
                        paste0(str_split(pheno_path[i],
                                         "/",
                                         simplify = TRUE)[,7],
                               ".",
                               str_split(pheno_path[i],
                                         "/",
                                         simplify = TRUE)[,8]))
  
  pheno <- full_join(pheno,
                     a,
                     by = "genome_id")
  
}

write_csv(pheno,
          "data/pheno_full.csv")
