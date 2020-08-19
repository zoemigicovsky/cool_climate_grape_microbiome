rm(list=ls(all.names=TRUE))

setwd("/home/gavin/github_repos/root_depth/scripts/")

source("root_depth_project_functions.R")

library(matrixStats)
library(openxlsx)

bacteria_meta <- read.table("../data/metadata/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)


bacteria_taxa <- read.table("../data/ASV_tables/bacteria/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

bacteria_taxa_breakdown$asv <- paste(bacteria_taxa_breakdown$species, ";ASV_", rownames(bacteria_taxa_breakdown), sep="")

bacteria_ASVs <- read.table("../data/ASV_tables/bacteria/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]



# Run same commands, but for fungi:
fungi_ASVs <- read.table("../data/ASV_tables/fungi/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]


fungi_taxa <- read.table("../data/ASV_tables/fungi/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

fungi_taxa_breakdown$asv <- paste(fungi_taxa_breakdown$species, ";ASV_", rownames(fungi_taxa_breakdown), sep="")

fungi_meta <- read.table("../data/metadata/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_meta <- fungi_meta[which(rownames(fungi_meta) %in% colnames(fungi_ASVs)), ]


bacteria_samples <- list()
bacteria_samples[["grape_root"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape"), ])
bacteria_samples[["grape_soil"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "soil" & bacteria_meta$species == "grape"), ])
bacteria_samples[["cover_crop_root"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "cover_crop"), ])

bacteria_samples[["grape_root_depth_0_15"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$root_depth == "depth_0_15"), ])
bacteria_samples[["grape_root_depth_15_30"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$root_depth == "depth_15_30"), ])
bacteria_samples[["grape_root_depth_30_50"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$root_depth == "depth_30_50"), ])

bacteria_samples[["grape_root_new_york_muscat"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$rootstock == "new_york_muscat"), ])
bacteria_samples[["grape_root_c3309"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$rootstock == "c3309"), ])
bacteria_samples[["grape_root_riparia_gloire"]] <- rownames(bacteria_meta[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$rootstock == "riparia_gloire"), ])




fungi_samples <- list()
fungi_samples[["grape_root"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape"), ])
fungi_samples[["grape_soil"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "soil" & fungi_meta$species == "grape"), ])
fungi_samples[["cover_crop_root"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "cover_crop"), ])

fungi_samples[["grape_root_depth_0_15"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$root_depth == "depth_0_15"), ])
fungi_samples[["grape_root_depth_15_30"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$root_depth == "depth_15_30"), ])
fungi_samples[["grape_root_depth_30_50"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$root_depth == "depth_30_50"), ])

fungi_samples[["grape_root_new_york_muscat"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$rootstock == "new_york_muscat"), ])
fungi_samples[["grape_root_c3309"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$rootstock == "c3309"), ])
fungi_samples[["grape_root_riparia_gloire"]] <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$rootstock == "riparia_gloire"), ])





# Get relative abundance tables sorted for each taxonomic level.
for(sample_group in names(bacteria_samples)) {

  tax_levels <- c("asv", "species", "genus", "family", "order", "class", "phylum")
  tax_levels <- rev(tax_levels)
  
  taxa_metrics <- list()
  
  for(kingdom in c("bacteria", "fungi")) {
    
    for(tax_level in tax_levels) {
      
      label <- paste(kingdom, tax_level, sep="_")
      
      if(kingdom == "bacteria") {
        samples <- bacteria_samples[[sample_group]]
        asv_table <- bacteria_ASVs[, samples]
        taxa_breakdown <- bacteria_taxa_breakdown
      } else {
        samples <- fungi_samples[[sample_group]]
        asv_table <- fungi_ASVs[, samples]
        taxa_breakdown <- fungi_taxa_breakdown
      }
      
      asv_table <- data.frame(sweep(asv_table, 2, colSums(asv_table), '/')) * 100
      
      empty_rows <- which(rowSums(asv_table) == 0)
      if(length(empty_rows) > 0) { asv_table <- asv_table[-empty_rows, ] }
      
      asv_table$tax_level <- taxa_breakdown[rownames(asv_table), tax_level]
      taxa_table <- aggregate(. ~ tax_level, data=asv_table, FUN=sum)
      rownames(taxa_table) <- taxa_table$tax_level
      
      taxa_table <- taxa_table[, -which(colnames(taxa_table) == "tax_level")]
      
      taxa_table <- taxa_table[order(rowMeans(taxa_table), decreasing=TRUE), ]
      
      taxa_metrics[[label]] <- data.frame(matrix(NA, nrow=nrow(taxa_table), ncol=5))
      rownames(taxa_metrics[[label]]) <- rownames(taxa_table)
      colnames(taxa_metrics[[label]]) <- c("taxon", "prevalence", "mean_RA_orig_order", "mean_RA", "CoefVar_RA")
      
      taxa_metrics[[label]]$taxon <- rownames(taxa_metrics[[label]])
      
      taxa_metrics[[label]]$mean_RA_orig_order <- 1:nrow(taxa_metrics[[label]])
      
      if(length(grep("Unclassified$", rownames(taxa_table))) > 0) {
        taxa_metrics[[label]] <- taxa_metrics[[label]][-grep("Unclassified$", rownames(taxa_table)), ]
        taxa_table <- taxa_table[-grep("Unclassified$", rownames(taxa_table)), ]
      }
      
      taxa_metrics[[label]]$prevalence <- (rowSums(taxa_table > 0) / ncol(taxa_table)) * 100
      taxa_metrics[[label]]$mean_RA <- rowMeans(taxa_table)
      taxa_metrics[[label]]$CoefVar_RA <- rowSds(as.matrix(taxa_table)) / taxa_metrics[[label]]$mean_RA
    }
  }
  
  outfile <- paste("../data/group_sorted_relabun/all_levels_", sample_group, "_metrics.xlsx", sep="")
  write.xlsx(taxa_metrics, file = outfile)

}



# Sanity checks on individual relative abundances to check above code.
# The expected mean relabun and CV were printed out to do quick checks by eye.
rm(list=ls(all.names=TRUE))

source("root_depth_project_functions.R")

library(matrixStats)

bacteria_meta <- read.table("../data/metadata/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)


bacteria_taxa <- read.table("../data/ASV_tables/bacteria/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

bacteria_ASVs <- read.table("../data/ASV_tables/bacteria/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]



fungi_ASVs <- read.table("../data/ASV_tables/fungi/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]


fungi_taxa <- read.table("../data/ASV_tables/fungi/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

fungi_taxa_breakdown$asv <- paste(fungi_taxa_breakdown$species, ";ASV_", rownames(fungi_taxa_breakdown), sep="")

fungi_meta <- read.table("../data/metadata/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_meta <- fungi_meta[which(rownames(fungi_meta) %in% colnames(fungi_ASVs)), ]



# test_asvs <- bacteria_ASVs
# test_samples <- rownames(bacteria_meta)[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape")]
# test_taxa <- bacteria_taxa_breakdown
# test_level <- "genus"

# test_asvs <- bacteria_ASVs
# test_samples <- rownames(bacteria_meta)[which(bacteria_meta$tissue == "root" & bacteria_meta$species == "grape" & bacteria_meta$root_depth == "depth_15_30")]
# test_taxa <- bacteria_taxa_breakdown
# test_level <- "phylum"

# test_asvs <- fungi_ASVs
# test_samples <- rownames(fungi_meta[which(fungi_meta$tissue == "root" & fungi_meta$species == "grape" & fungi_meta$rootstock == "riparia_gloire"), ])
# test_taxa <- fungi_taxa_breakdown
# test_level <- "family"


test_asvs <- bacteria_ASVs
test_samples <- rownames(bacteria_meta)[which(bacteria_meta$tissue == "soil" & bacteria_meta$species == "grape")]
test_taxa <- bacteria_taxa_breakdown
test_level <- "class"

test_asvs <- test_asvs[, test_samples]

test_asvs <- data.frame(sweep(test_asvs, 2, colSums(test_asvs), '/')) * 100

test_asvs$taxon <- test_taxa[rownames(test_asvs), test_level]

test_aggregated <- aggregate(. ~ taxon, data=test_asvs, FUN = sum)

rownames(test_aggregated) <- test_aggregated$taxon
test_aggregated <- test_aggregated[, -which(colnames(test_aggregated) == "taxon")]


test_aggregated <- data.frame(sweep(test_aggregated, 2, colSums(test_aggregated), '/')) * 100

sorted_order <- names(sort(rowSums(test_aggregated), decreasing = TRUE))

test_aggregated_sorted <- test_aggregated[sorted_order, ]

rowMeans(test_aggregated_sorted)[1:10]

(rowSds(as.matrix(test_aggregated_sorted)) / rowMeans(test_aggregated_sorted))[1:10]
