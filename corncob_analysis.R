library(tidyverse)
library(corncob)
library(phyloseq)

###BACTERIA####
bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
#Create column for the plant id 
bacteria_meta$sample_id <- rownames(bacteria_meta)

bacteria_meta <- bacteria_meta %>% mutate(other_info=str_replace(other_info, "NYM Own Root ", ""),other_info=str_replace(other_info, "NYM 3309", ""),other_info=str_replace(other_info, "NYM Rg", ""),other_info=str_replace(other_info, "0-15cm", ""),other_info=str_replace(other_info, "15-30cm", ""),other_info=str_replace(other_info, "30-50cm", ""),other_info=str_replace(other_info, "NYM", ""),other_info=str_trim(other_info))

bacteria_meta <- as.data.frame(bacteria_meta)
rownames(bacteria_meta) <- bacteria_meta$sample_id
bacteria_meta <- bacteria_meta[,-ncol(bacteria_meta)]

#load in grape root only ASVs
bacteria_ASVs <- read.table("bacteria/dada2_output_exported_grape/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

#Keep only meta data for grape roots
bacteria_meta <- bacteria_meta[colnames(bacteria_ASVs),]

table(colnames(bacteria_ASVs)== rownames(bacteria_meta))

#use the ASV table and metadata table to create phyloseq otu table and sample data, then combine into a phyloseq object
otu <- otu_table(bacteria_ASVs, taxa_are_rows = TRUE)
sam_data <- sample_data(bacteria_meta, errorIfNULL = TRUE)
phylo <-merge_phyloseq(otu,sam_data)

set.seed(38)

#https://rdrr.io/a/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd

da_analysis <- differentialTest(formula = ~rootstock+root_depth,
                 phi.formula = ~1,
                 formula_null = ~rootstock+ root_depth,
                 phi.formula_null = ~1,
                 test= "Wald", boot = FALSE,
                 data = phylo,
                 fdr_cutoff = 0.05)

da_analysis$significant_taxa

plot(da_analysis2)
