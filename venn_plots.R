# Read in files and determine genera overlapping in cover, soil, and root for bacteria and fungal datasets.

rm(list=ls(all=TRUE))

library("VennDiagram")

# Note that these two lines are specific to Vulcan:
setwd("/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/")
source("/home/gavin/github_repos/root_depth/root_depth_project_functions.R")

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

bacteria_meta$group <- bacteria_meta$tissue
bacteria_meta$group[which(bacteria_meta$species == "cover_crop")] <- "root (cover)"

bacteria_ASVs <- read.table("bacteria/dada2_output_exported/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

bacteria_taxa <- read.table("bacteria/taxa/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

# Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
bacteria_root_grape_samples <- rownames(bacteria_meta)[which(bacteria_meta$group == "root")]
bacteria_root_cover_samples <- rownames(bacteria_meta)[which(bacteria_meta$group == "root (cover)")]
bacteria_soil_grape_samples <- rownames(bacteria_meta)[which(bacteria_meta$group == "soil")]

root_grape_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_root_grape_samples]) > 0)]
root_grape_genera <- unique(bacteria_taxa_breakdown[root_grape_ASVs, "genus"])

root_cover_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_root_cover_samples]) > 0)]
root_cover_genera <- unique(bacteria_taxa_breakdown[root_cover_ASVs, "genus"])

soil_grape_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_soil_grape_samples]) > 0)]
soil_grape_genera <- unique(bacteria_taxa_breakdown[soil_grape_ASVs, "genus"])

bacteria_root_grape_genus_count <- length(root_grape_genera)
bacteria_root_cover_genus_count <- length(root_cover_genera)
bacteria_soil_grape_genus_count <- length(soil_grape_genera)
bacteria_root_grape_cover_genus_count <- length(which(root_grape_genera %in% root_cover_genera))
bacteria_root_cover_soil_grape_genus_count <- length(which(root_cover_genera %in% soil_grape_genera))
bacteria_root_soil_grape_genus_count <- length(which(root_grape_genera %in% soil_grape_genera))
bacteria_all3_genus_count_TMP <- root_grape_genera[which(root_grape_genera %in% soil_grape_genera)]
bacteria_all3_genus_count <- length(bacteria_all3_genus_count_TMP[which(bacteria_all3_genus_count_TMP %in% root_cover_genera)])

# Saved as 7x7 inches
draw.triple.venn(area1=bacteria_root_grape_genus_count,
                 area2=bacteria_root_cover_genus_count,
                 area3=bacteria_soil_grape_genus_count,
                 n12 = bacteria_root_grape_cover_genus_count,
                 n23 = bacteria_root_cover_soil_grape_genus_count,
                 n13 = bacteria_root_soil_grape_genus_count,
                 n123 = bacteria_all3_genus_count,
                 category = c("Grape roots", "Cover crop roots", "Grape soil"),
                 scaled=TRUE,
                 fill = c("blue", "red", "yellow"))
                 

#### Run the same commands, but for fungi:
rm(list=ls(all=TRUE))

library("VennDiagram")

# Note that these two lines are specific to Vulcan:
setwd("/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/")
source("/home/gavin/github_repos/root_depth/root_depth_project_functions.R")

fungi_ASVs <- read.table("fungi/dada2_output_exported/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

fungi_taxa <- read.table("fungi/taxa/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_meta$group <- fungi_meta$tissue
fungi_meta$group[which(fungi_meta$species == "cover_crop")] <- "root (cover)"

# Need to subset metadata to only samples in BIOM table.
fungi_meta <- fungi_meta[which(rownames(fungi_meta) %in% colnames(fungi_ASVs)), ]

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

# Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
fungi_root_grape_samples <- rownames(fungi_meta)[which(fungi_meta$group == "root")]
fungi_root_cover_samples <- rownames(fungi_meta)[which(fungi_meta$group == "root (cover)")]
fungi_soil_grape_samples <- rownames(fungi_meta)[which(fungi_meta$group == "soil")]

root_grape_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_root_grape_samples]) > 0)]
root_grape_genera <- unique(fungi_taxa_breakdown[root_grape_ASVs, "genus"])

root_cover_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_root_cover_samples]) > 0)]
root_cover_genera <- unique(fungi_taxa_breakdown[root_cover_ASVs, "genus"])

soil_grape_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_soil_grape_samples]) > 0)]
soil_grape_genera <- unique(fungi_taxa_breakdown[soil_grape_ASVs, "genus"])

fungi_root_grape_genus_count <- length(root_grape_genera)
fungi_root_cover_genus_count <- length(root_cover_genera)
fungi_soil_grape_genus_count <- length(soil_grape_genera)
fungi_root_grape_cover_genus_count <- length(which(root_grape_genera %in% root_cover_genera))
fungi_root_cover_soil_grape_genus_count <- length(which(root_cover_genera %in% soil_grape_genera))
fungi_root_soil_grape_genus_count <- length(which(root_grape_genera %in% soil_grape_genera))
fungi_all3_genus_count_TMP <- root_grape_genera[which(root_grape_genera %in% soil_grape_genera)]
fungi_all3_genus_count <- length(fungi_all3_genus_count_TMP[which(fungi_all3_genus_count_TMP %in% root_cover_genera)])

# Saved as 7x7 inches
draw.triple.venn(area1=fungi_root_grape_genus_count,
                 area2=fungi_root_cover_genus_count,
                 area3=fungi_soil_grape_genus_count,
                 n12 = fungi_root_grape_cover_genus_count,
                 n23 = fungi_root_cover_soil_grape_genus_count,
                 n13 = fungi_root_soil_grape_genus_count,
                 n123 = fungi_all3_genus_count,
                 category = c("Grape roots", "Cover crop roots", "Grape soil"),
                 scaled=TRUE,
                 fill = c("blue", "red", "yellow"))

