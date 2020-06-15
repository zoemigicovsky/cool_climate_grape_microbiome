### For several tables the ordering of taxa from greatest to lowest relative abundance is of interest. 
### This is specifically of interest in the grape root samples.
### Goal here was to aggregate the rows by the different taxa levels and to output these tables.
### Also output a table with the taxa ids in sorted order for each column and for bacteria/fungi separately.
### Note that these analyses were based on the previously rarified biom tables.

rm(list=ls(all.names=TRUE))

library(reshape2)
library(tidyverse)

# Change this line for different systems:
path2repo <- "/home/gavin/github_repos/root_depth"

setwd(path2repo)
source("scripts/root_depth_project_functions.R")

bacteria_meta <- read.table("data/metadata/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

bacteria_graperoot_samples <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape"), ])


bacteria_ASVs <- read.table("data/ASV_tables/bacteria/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs_graperoot <- bacteria_ASVs[, bacteria_graperoot_samples]

bacteria_ASVs_graperoot_rel <- data.frame(sweep(bacteria_ASVs_graperoot, 2, colSums(bacteria_ASVs_graperoot), '/')) * 100

bacteria_ASVs_graperoot_rel <- bacteria_ASVs_graperoot_rel[-which(rowSums(bacteria_ASVs_graperoot_rel) == 0), ]

bacteria_taxa <- read.table("data/ASV_tables/bacteria/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

bacteria_taxa_breakdown$asv <- paste(bacteria_taxa_breakdown$species, ";ASV_", rownames(bacteria_taxa_breakdown), sep="")



# Run same commands, but for fungi:
fungi_ASVs <- read.table("data/ASV_tables/fungi/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

#Convert to relative abundance (i.e. sum to 100%) - the "sweep" function is helpful for this.
fungi_asv_abun <- data.frame(sweep(fungi_ASVs, 2, colSums(fungi_ASVs), '/')) * 100


fungi_taxa <- read.table("data/ASV_tables/fungi/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

fungi_meta <- read.table("data/metadata/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_graperoot_samples <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape"), ])

fungi_graperoot_samples <- fungi_graperoot_samples[which(fungi_graperoot_samples %in% colnames(fungi_ASVs))]

fungi_ASVs_graperoot <- fungi_ASVs[, fungi_graperoot_samples]

fungi_ASVs_graperoot_rel <- fungi_ASVs_graperoot

fungi_ASVs_graperoot_rel <- data.frame(sweep(fungi_ASVs_graperoot, 2, colSums(fungi_ASVs_graperoot), '/')) * 100

fungi_ASVs_graperoot_rel <- fungi_ASVs_graperoot_rel[-which(rowSums(fungi_ASVs_graperoot_rel) == 0), ]

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

fungi_taxa_breakdown$asv <- paste(fungi_taxa_breakdown$species, ";ASV_", rownames(fungi_taxa_breakdown), sep="")



# Get relative abundance tables sorted for each taxonomic level.
tax_levels <- c("asv", "species", "genus", "family", "order", "class", "phylum")
tax_levels <- rev(tax_levels)

graperoot_sorted_taxa <- list()

graperoot_sorted_collapsed_relabun <- list()

for(tax_level in tax_levels) {
  
  bacteria_label <- paste("bacteria", tax_level, sep="_")
  fungi_label <- paste("fungi", tax_level, sep="_")
  
  tmp_bacteria <- bacteria_ASVs_graperoot_rel
  tmp_bacteria$tax_level <- bacteria_taxa_breakdown[rownames(tmp_bacteria), tax_level]
  tmp_bacteria <- aggregate(. ~ tax_level, data=tmp_bacteria, FUN=sum)
  tmp_bacteria <- tmp_bacteria[order(rowSums(tmp_bacteria[, 2:ncol(tmp_bacteria)]), decreasing=TRUE), ]
  rownames(tmp_bacteria) <- tmp_bacteria$tax_level
  tmp_bacteria <- tmp_bacteria[, -1]
  
  graperoot_sorted_collapsed_relabun[[bacteria_label]] <- tmp_bacteria
  graperoot_sorted_taxa[[bacteria_label]] <- rownames(tmp_bacteria)
  
  tmp_fungi <- fungi_ASVs_graperoot_rel
  tmp_fungi$tax_level <- fungi_taxa_breakdown[rownames(tmp_fungi), tax_level]
  tmp_fungi <- aggregate(. ~ tax_level, data=tmp_fungi, FUN=sum)
  tmp_fungi <- tmp_fungi[order(rowSums(tmp_fungi[, 2:ncol(tmp_fungi)]), decreasing=TRUE), ]
  rownames(tmp_fungi) <- tmp_fungi$tax_level
  tmp_fungi <- tmp_fungi[, -1]
  
  graperoot_sorted_collapsed_relabun[[fungi_label]] <- tmp_fungi
  graperoot_sorted_taxa[[fungi_label]] <- rownames(tmp_fungi)
  
  
}

# Write output tables
sorted_taxa_lists <- data.frame(matrix(NA, nrow=length(graperoot_sorted_taxa$bacteria_asv), ncol=14))
colnames(sorted_taxa_lists) <- c(paste("bacteria", tax_levels, sep="_"), paste("fungi", tax_levels, sep="_"))

for(taxa_set in colnames(sorted_taxa_lists)) {
 
  taxa_order <- graperoot_sorted_taxa[[taxa_set]]
  
  if(length(taxa_order) < nrow(sorted_taxa_lists)) {
    taxa_order <- c(taxa_order, rep(NA, nrow(sorted_taxa_lists) - length(taxa_order)))
  }
   
  sorted_taxa_lists[, taxa_set] <- taxa_order
  
  write.table(x=graperoot_sorted_collapsed_relabun[[taxa_set]],
              file=paste("data/grape_root_sorted_relabun/", taxa_set, "_sorted_relabun.tsv", sep=""),
              row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")
  
}

write.table(x = sorted_taxa_lists,
            file = "data/grape_root_sorted_relabun/taxa_sorted_by_relabun.tsv",
            row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")


# Sanity checks:
# tmp <- graperoot_sorted_collapsed_relabun$bacteria_asv
# rownames(tmp) <- gsub("^.*ASV_", "", rownames(tmp))
# all.equal(tmp, bacteria_ASVs_graperoot_rel[rownames(tmp), ])
# 
# tmp <- graperoot_sorted_collapsed_relabun$fungi_asv
# rownames(tmp) <- gsub("^.*ASV_", "", rownames(tmp))
# all.equal(tmp, fungi_ASVs_graperoot_rel[rownames(tmp), ])
# 
# 
# test_asvs <- rownames(fungi_taxa_breakdown[which(fungi_taxa_breakdown$order == "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales"), ])
# test_asvs <- test_asvs[which(test_asvs %in% rownames(fungi_ASVs_graperoot_rel))]
# colSums(fungi_ASVs_graperoot_rel[test_asvs, ])

