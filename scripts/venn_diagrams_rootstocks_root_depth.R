# Make venn diagrams to show overlapping genera for rootstocks and root depths (based on bacteria and fungi). 
# Note that RDS files of rarefied tables were made in the figure1 Rscript.

rm(list=ls(all=TRUE))

library(extrafont)
library(cowplot)
library(VennDiagram)
source("scripts/root_depth_project_functions.R")

###BACTERIA ROOTSTOCK####

bacteria_meta <- read.table("data/metadata/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs (subsampled in advance).
bacteria_ASVs_1176subsample_RDS <- paste("data/intermediate_RDS/bacteria_ASVs_1176subsample.rds")
bacteria_ASVs_1176subsample <- readRDS(bacteria_ASVs_1176subsample_RDS)

# Only keep grape root samples.
bacteria_meta <- bacteria_meta[which(bacteria_meta$species == "grape" & bacteria_meta$tissue == "root"), ]
bacteria_ASVs_1176subsample <- bacteria_ASVs_1176subsample[ , rownames(bacteria_meta)]

#Get taxonomy information
bacteria_taxa <- read.table("data/ASV_tables/bacteria/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

bacteria_rootstock_venn <- threeWayVennPercentGenus(metadata=bacteria_meta,
                                                meta_col="rootstock",
                                                asv_abun=bacteria_ASVs_1176subsample,
                                                taxa_df=bacteria_taxa_breakdown,
                                                meta_cat=c("new_york_muscat", "c3309", "riparia_gloire"),
                                                labels=c("Ungrafted", "3309 C", "Riparia Gloire"),
                                                colours=c("#009E73", "#E69F00", "#56B4E9"))

###FUNGI ROOTSTOCK####

fungi_meta <- read.table("data/metadata/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs (subsampled in advance)
fungi_ASVs_1176subsample_RDS <- paste(path2repo, "data/intermediate_RDS/fungi_ASVs_1176subsample.rds", sep="/")
fungi_ASVs_1176subsample <- readRDS("data/intermediate_RDS/fungi_ASVs_1176subsample.rds")


# Only keep grape root samples.
fungi_meta <- fungi_meta[which(fungi_meta$species == "grape" & fungi_meta$tissue == "root"), ]
fungi_overlapping_samples <- colnames(fungi_ASVs_1176subsample)[which(colnames(fungi_ASVs_1176subsample) %in% rownames(fungi_meta))]
fungi_ASVs_1176subsample <- fungi_ASVs_1176subsample[, fungi_overlapping_samples]
fungi_meta <- fungi_meta[fungi_overlapping_samples, ]

#Get taxonoomy information
fungi_taxa <- read.table("data/ASV_tables/fungi/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- qiime2_taxa_breakdown(fungi_taxa)

fungi_rootstock_venn <- threeWayVennPercentGenus(metadata=fungi_meta,
                                             meta_col="rootstock",
                                             asv_abun=fungi_ASVs_1176subsample,
                                             taxa_df=fungi_taxa_breakdown,
                                             meta_cat=c("new_york_muscat", "c3309", "riparia_gloire"),
                                             labels=c("Ungrafted", "3309 C", "Riparia Gloire"),
                                             colours=c("#009E73", "#E69F00", "#56B4E9"))

pdf("figures/rootstock_venn_diagram.pdf", width=6.5, height=3, family="Arial")
plot_grid(grobTree(bacteria_rootstock_venn), grobTree(fungi_rootstock_venn), labels="AUTO", label_size=25)
dev.off()

#Make venn diagrams to show overlapping genera for root_depths and root depths (based on bacteria and fungi). 

###BACTERIA ROOT DEPTH####

bacteria_meta <- read.table("data/metadata/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs (subsampled in advance).
bacteria_ASVs_1176subsample_RDS <- paste("data/intermediate_RDS/bacteria_ASVs_1176subsample.rds")
bacteria_ASVs_1176subsample <- readRDS(bacteria_ASVs_1176subsample_RDS)

# Only keep grape root samples.
bacteria_meta <- bacteria_meta[which(bacteria_meta$species == "grape" & bacteria_meta$tissue == "root"), ]
bacteria_ASVs_1176subsample <- bacteria_ASVs_1176subsample[ , rownames(bacteria_meta)]

#Get taxonomy information
bacteria_taxa <- read.table("data/ASV_tables/bacteria/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)


bacteria_root_depth_venn <- threeWayVennPercentGenus(metadata=bacteria_meta,
                                                meta_col="root_depth",
                                                asv_abun=bacteria_ASVs_1176subsample,
                                                taxa_df=bacteria_taxa_breakdown,
                                                meta_cat=c("depth_0_15", "depth_15_30", "depth_30_50"),
                                                labels=c("0-15 cm", "15-30 cm", "30-50 cm"),
                                                colours=c("#009E73", "#E69F00", "#56B4E9"))

###FUNGI ROOT DEPTH####

fungi_meta <- read.table("data/metadata/root_depth_fungi_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs (subsampled above).
fungi_ASVs_1176subsample <- readRDS(fungi_ASVs_1176subsample_RDS)

# Only keep grape root samples.
fungi_meta <- fungi_meta[which(fungi_meta$species == "grape" & fungi_meta$tissue == "root"), ]
fungi_overlapping_samples <- colnames(fungi_ASVs_1176subsample)[which(colnames(fungi_ASVs_1176subsample) %in% rownames(fungi_meta))]
fungi_ASVs_1176subsample <- fungi_ASVs_1176subsample[, fungi_overlapping_samples]
fungi_meta <- fungi_meta[fungi_overlapping_samples, ]

#Get taxonoomy information
fungi_taxa <- read.table("data/ASV_tables/fungi/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- qiime2_taxa_breakdown(fungi_taxa)

fungi_root_depth_venn <- threeWayVennPercentGenus(metadata=fungi_meta,
                                             meta_col="root_depth",
                                             asv_abun=fungi_ASVs_1176subsample,
                                             taxa_df=fungi_taxa_breakdown,
                                             meta_cat=c("depth_0_15", "depth_15_30", "depth_30_50"),
                                             labels=c("0-15 cm", "15-30 cm", "30-50 cm"),
                                             colours=c("#009E73", "#E69F00", "#56B4E9"))

pdf("figures/root_depth_venn_diagram.pdf", width=6.5, height=3, family="Arial")
plot_grid(grobTree(bacteria_root_depth_venn), grobTree(fungi_root_depth_venn), labels="AUTO", label_size=25)
dev.off()
