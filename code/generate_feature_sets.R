#LIBRARIES ----
library(tidyverse)

#FEATURES ----
cog_subset <- read_csv("data/cog_subset.csv",
                       col_names = TRUE)
kegg_subset <- read_csv("data/kegg_subset.csv",
                        col_names = TRUE)
key_terms <- read_csv("data/key_terms.csv",
                      col_names = TRUE)
gene_terms <- read_csv("data/gene_terms.csv",
                       col_names = TRUE)

annots_full <- read_csv("../../../data/complete_eggnog_annots.csv") %>%
  mutate(locus_tag = str_extract(seed_eggNOG_ortholog, "CD630_....."))

ncbi_gene_files <- list.files("../../data/reference_gene_data",
                              pattern = ".csv",
                              full.names = TRUE)

ncbi_gene_pull <- data.table::rbindlist(lapply(ncbi_gene_files,
                                               read_csv))

ncbi_gene <- ncbi_gene_pull %>%
  mutate(`Locus tag` = gsub("RS", "", `Locus tag`))

colnames(ncbi_gene) <- c("Genome_loc",
                         colnames(ncbi_gene)[2:6],
                         "Gene_name",
                         "locus_tag",
                         "protein_product",
                         "length",
                         "description")

combo_annots_full <- full_join(annots_full,
                               ncbi_gene,
                               by = intersect(colnames(annots_full),
                                              colnames(ncbi_gene))) %>%
  distinct()

#SEARCH DESCRIPTIONS ----
desc_index <- lapply(key_terms$Key_Terms,
                     function(x){
                       
                       desc1 <- grep(x, combo_annots_full$`eggNOG free text desc.`)
                       desc2 <- grep(x, combo_annots_full$description)
                       
                       out <- unique(c(desc1, desc2))
                       
                       return(out)
                       
                     })
names(desc_index) <- key_terms$Key_Terms
desc_index_unique <- unique(unlist(desc_index))

#SEARCH GENE NAMES ----
gene_index <- lapply(gene_terms$Gene_Names,
                     function(x){
                       
                       gene1 <- grep(x, combo_annots_full$Preferred_name)
                       gene2 <- grep(x, combo_annots_full$Gene_name)
                       
                       out <- unique(c(gene1, gene2))
                       
                       return(out)
                       
                     })
names(gene_index) <- gene_terms$Gene_Names
gene_index_unique <- unique(unlist(gene_index))

#SEARCH COGS ----
cog_index <- lapply(cog_subset$COG_Letter,
                    function(x){
                      
                      out <- grep(x, combo_annots_full$`COG Functional cat.`)
                      
                      return(out)
                    })
names(cog_index) <- cog_subset$COG_Letter
cog_index_unique <- unique(unlist(cog_index))

#SEARCH KEGG ----
kegg_index <- lapply(kegg_subset$Map_ID,
                     function(x){
                       
                       out <- grep(x, combo_annots_full$KEGG_Pathway)
                       
                       return(out)
                     })
names(kegg_index) <- kegg_subset$Map_ID
kegg_index_unique <- unique(unlist(kegg_index))

#COMBINED UNIQUE INDICES ----
index <- unique(c(desc_index_unique,
                  gene_index_unique,
                  cog_index_unique,
                  kegg_index_unique))

#CLEAN ANNOTS ----
combo_annots_sub <- combo_annots_full[index,]

panaroo_annots <- combo_annots_sub %>%
  filter(dataset == "panaroo")

core_annots <- combo_annots_sub %>%
  filter(grepl("CD630_", locus_tag))

#GENOS ----
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
