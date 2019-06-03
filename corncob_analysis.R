library(tidyverse)
library(corncob)
library(phyloseq)

###BACTERIA####
bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata_grape_vine.txt",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

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

da_analysis_rootstock <- differentialTest(formula = ~rootstock+root_depth,
                                phi.formula = ~rootstock+root_depth,
                                formula_null = ~root_depth,
                                phi.formula_null = ~rootstock+root_depth,
                                test= "Wald", boot = FALSE,
                                data = phylo,
                                fdr_cutoff = 0.05)

da_analysis_rootstock$significant_taxa

plot(da_analysis_rootstock)


da_analysis_root_depth <- differentialTest(formula = ~rootstock+root_depth,
                                           phi.formula = ~rootstock+root_depth,
                                           formula_null = ~rootstock,
                                           phi.formula_null = ~rootstock+root_depth,
                                           test= "Wald", boot = FALSE,
                                           data = phylo,
                                           fdr_cutoff = 0.05)

da_analysis_root_depth$significant_taxa

plot(da_analysis_root_depth)