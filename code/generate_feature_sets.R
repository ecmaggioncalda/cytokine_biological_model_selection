#LIBRARIES ----
library(tidyverse)

#DATA ----
load("../../2023_02_14_biological_model_explore/data/eggNOG_kegg_cog_annots_2023_04_11.RData")

#CLEAN ANNOTS ----
index_annots_count <- index_df %>%
  select(`#query_name`,
         seed_eggNOG_ortholog,
         locus_tag,
         contains("hits")) %>%
  distinct() %>%
  group_by(seed_eggNOG_ortholog,
           locus_tag) %>%
  summarize(across(where(is.numeric), ~sum(.x)))

length(unique(index_annots_count$locus_tag[!is.na(index_annots_count$locus_tag)]))

panaroo_annots <- index_df %>%
  filter(dataset == "panaroo") %>%
  select(`#query_name`) %>%
  drop_na() %>%
  distinct()

length(unique(panaroo_annots$`#query_name`))

core_annots <- index_df %>%
  select(locus_tag) %>%
  drop_na() %>%
  distinct()

length(unique(core_annots$locus_tag))

#PAN ----
pan_df_full <- read_delim("../../data/pan_mat.tsv")

pan_variants_full <- gsub("$", "_", pan_df_full$variant)

pan_variants_index <- sapply(pan_variants_full,
                             function(x){
                               
                               any(x == panaroo_annots$`#query_name`)
                               
                             })

table(pan_variants_index)

pan_df_sub <- pan_df_full[pan_variants_index,]

#STRUCT ----
struct_df_full <- read_delim("../../data/pan_struct_mat.tsv")

struct_names_pull <- struct_df_full %>%
  select(variant) %>%
  mutate(variant = gsub("\\.\\.\\.", "~", variant)) %>%
  deframe()

struct_names_split <- str_split(struct_names_pull, "\\.", simplify = TRUE) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~gsub("~", "\\.\\.\\.", .x)),
         across(everything(), ~gsub("$", "_", .x)))

struct_variants_index1 <- sapply(struct_names_split[,1],
                                 function(x){
                                   
                                   any(x == panaroo_annots$`#query_name`)
                                   
                                 })

struct_variants_index2 <- sapply(struct_names_split[,2],
                                 function(x){
                                   
                                   any(x == panaroo_annots$`#query_name`)
                                   
                                 })

struct_variants_index3 <- sapply(struct_names_split[,3],
                                 function(x){
                                   
                                   any(x == panaroo_annots$`#query_name`)
                                   
                                 })

struct_variants_index <- rowSums(cbind(struct_variants_index1,
                                       struct_variants_index2,
                                       struct_variants_index3)) >= 1

table(struct_variants_index)

struct_df_sub <- struct_df_full[struct_variants_index,]

#CORE ----
core_df_full <- read_delim("../../data/core_mat_sift.tsv")
core_df_annots_full <- read_csv("../../data/annots_mat_sift.csv")

core_variants_index <- sapply(core_df_annots_full$locus_tag,
                              function(x){
                                
                                any(x == core_annots$locus_tag)
                                
                              })

table(core_variants_index)

core_df_sub <- core_df_full[core_variants_index,]

#GENE ----
gene_df_full <- read_delim("../../data/gene_mat.tsv")

gene_variants_index <- sapply(gene_df_full$variant,
                              function(x){
                                
                                any(x == core_annots$locus_tag)
                                
                              })

table(gene_variants_index)

gene_df_sub <- gene_df_full[gene_variants_index,]

#COMBINED ----
geno_df_sub <- bind_rows(core_df_sub,
                         gene_df_sub,
                         pan_df_sub,
                         struct_df_sub)
dim(geno_df_sub)

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
