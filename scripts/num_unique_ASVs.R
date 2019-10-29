### Script to explore whether there is evidence of more rare ASVs in fungi root samples compared to bacterial root samples.
### There was no evidence for this once there were subsampled to the same read depth.

rm(list=ls(all=TRUE))

library(ggplot2)
library(vegan)

# Change this line for different systems:
path2repo <- "/home/gavin/github_repos/root_depth"

setwd(path2repo)
source("scripts/root_depth_project_functions.R")

root_meta <- read.table("data/metadata/root_depth_bacteria_metadata.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
root_meta$group <- root_meta$tissue
root_meta$group[which(root_meta$species == "cover_crop")] <- "root (cover)"

bacteria_ASVs_1176subsample_RDS <- paste("data/intermediate_RDS/bacteria_ASVs_1176subsample.rds")
bacteria_ASVs_1176subsample <- readRDS(bacteria_ASVs_1176subsample_RDS)

fungi_ASVs_1176subsample_RDS <- paste("data/intermediate_RDS/fungi_ASVs_1176subsample.rds")
fungi_ASVs_1176subsample <- readRDS(fungi_ASVs_1176subsample_RDS)

root_samples <- rownames(root_meta)[which(root_meta$group == "root")]
root_samples <- root_samples[which(root_samples %in% colnames(fungi_ASVs_1176subsample))]
root_samples <- root_samples[which(root_samples %in% colnames(bacteria_ASVs_1176subsample))]


bacteria_ASVs_1176subsample_root <- bacteria_ASVs_1176subsample[, root_samples]
fungi_ASVs_1176subsample_root <- fungi_ASVs_1176subsample[, root_samples]

bacteria_ASVs_1176subsample_root_rowSums <- rowSums(bacteria_ASVs_1176subsample_root > 0)
fungi_ASVs_1176subsample_root_rowSums <- rowSums(fungi_ASVs_1176subsample_root > 0)

bacteria_ASVs_1176subsample_root_rowSums <- bacteria_ASVs_1176subsample_root_rowSums[-which(bacteria_ASVs_1176subsample_root_rowSums == 0)]
fungi_ASVs_1176subsample_root_rowSums <- fungi_ASVs_1176subsample_root_rowSums[-which(fungi_ASVs_1176subsample_root_rowSums == 0)]

boxplot(bacteria_ASVs_1176subsample_root_rowSums,
        fungi_ASVs_1176subsample_root_rowSums,
        names = c("Bacteria", "Fungi"),
        ylab="Number of samples positive for ASV", col="grey")

wilcox.test(fungi_ASVs_1176subsample_root_rowSums, bacteria_ASVs_1176subsample_root_rowSums)

# Wilcoxon rank sum test with continuity correction
# data:  fungi_ASVs_1176subsample_root_rowSums and bacteria_ASVs_1176subsample_root_rowSums
# W = 172350, p-value = 0.2132
# alternative hypothesis: true location shift is not equal to 0
