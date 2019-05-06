#Make venn diagrams to show overlapping genera for rootstocks and root depths (based on bacteria and fungi). 

library("VennDiagram")

source("root_depth_project_functions.R")

###BACTERIA ROOTSTOCK####

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs
bacteria_ASVs <- read.table("bacteria/dada2_output_exported_grape/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

#Keep only meta data for grape roots
bacteria_meta <- bacteria_meta[colnames(bacteria_ASVs),]

table(colnames(bacteria_ASVs)== rownames(bacteria_meta))
#they match

#Get taxonoomy information
bacteria_taxa <- read.table("bacteria/taxa/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

# Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
bacteria_new_york_samples <- rownames(bacteria_meta)[which(bacteria_meta$rootstock == "new_york_muscat")]
bacteria_c3309_samples <- rownames(bacteria_meta)[which(bacteria_meta$rootstock == "c3309")]
bacteria_riparia_gloire_samples <- rownames(bacteria_meta)[which(bacteria_meta$rootstock == "riparia_gloire")]

new_york_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_new_york_samples]) > 0)]
new_york_genera <- unique(bacteria_taxa_breakdown[new_york_ASVs, "genus"])

c3309_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_c3309_samples]) > 0)]
c3309_genera <- unique(bacteria_taxa_breakdown[c3309_ASVs, "genus"])

rg_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_riparia_gloire_samples]) > 0)]
rg_genera <- unique(bacteria_taxa_breakdown[rg_ASVs, "genus"])

bacteria_new_york_genus_count <- length(new_york_genera)
bacteria_c3309_genus_count <- length(c3309_genera)
bacteria_rg_genus_count <- length(rg_genera)
bacteria_new_york_c3309_genus_count <- length(which(new_york_genera %in% c3309_genera))
bacteria_c3309_rg_genus_count <- length(which(c3309_genera %in% rg_genera))
bacteria_new_york_rg_genus_count <- length(which(new_york_genera %in% rg_genera))
bacteria_all3_genus_count_TMP <- new_york_genera[which(new_york_genera %in% rg_genera)]
bacteria_all3_genus_count <- length(bacteria_all3_genus_count_TMP[which(bacteria_all3_genus_count_TMP %in% c3309_genera)])

# Saved as 7x7 inches
bacteria_rootstock_venn <- draw.triple.venn(area1=bacteria_new_york_genus_count,
                                  area2=bacteria_c3309_genus_count,
                                  area3=bacteria_rg_genus_count,
                                  n12 = bacteria_new_york_c3309_genus_count,
                                  n23 = bacteria_c3309_rg_genus_count,
                                  n13 = bacteria_new_york_rg_genus_count,
                                  n123 = bacteria_all3_genus_count,
                                  category = c("new_york_muscat", "c3309", "riparia_gloire"),
                                  scaled=TRUE,
                                  fill = c( "#009E73","#E69F00","#56B4E9"))

###FUNGI ROOTSTOCK####

fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs
fungi_ASVs <- read.table("fungi/dada2_output_exported_grape/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

#Keep only meta data for grape roots
fungi_meta <- fungi_meta[colnames(fungi_ASVs),]

table(colnames(fungi_ASVs)== rownames(fungi_meta))
#they match

#Get taxonoomy information
fungi_taxa <- read.table("fungi/taxa/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- qiime2_taxa_breakdown(fungi_taxa)

# Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
fungi_new_york_samples <- rownames(fungi_meta)[which(fungi_meta$rootstock == "new_york_muscat")]
fungi_c3309_samples <- rownames(fungi_meta)[which(fungi_meta$rootstock == "c3309")]
fungi_riparia_gloire_samples <- rownames(fungi_meta)[which(fungi_meta$rootstock == "riparia_gloire")]

new_york_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_new_york_samples]) > 0)]
new_york_genera <- unique(fungi_taxa_breakdown[new_york_ASVs, "genus"])

c3309_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_c3309_samples]) > 0)]
c3309_genera <- unique(fungi_taxa_breakdown[c3309_ASVs, "genus"])

rg_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_riparia_gloire_samples]) > 0)]
rg_genera <- unique(fungi_taxa_breakdown[rg_ASVs, "genus"])

fungi_new_york_genus_count <- length(new_york_genera)
fungi_c3309_genus_count <- length(c3309_genera)
fungi_rg_genus_count <- length(rg_genera)
fungi_new_york_c3309_genus_count <- length(which(new_york_genera %in% c3309_genera))
fungi_c3309_rg_genus_count <- length(which(c3309_genera %in% rg_genera))
fungi_new_york_rg_genus_count <- length(which(new_york_genera %in% rg_genera))
fungi_all3_genus_count_TMP <- new_york_genera[which(new_york_genera %in% rg_genera)]
fungi_all3_genus_count <- length(fungi_all3_genus_count_TMP[which(fungi_all3_genus_count_TMP %in% c3309_genera)])

# Saved as 7x7 inches
fungi_rootstock_venn <- draw.triple.venn(area1=fungi_new_york_genus_count,
                 area2=fungi_c3309_genus_count,
                 area3=fungi_rg_genus_count,
                 n12 = fungi_new_york_c3309_genus_count,
                 n23 = fungi_c3309_rg_genus_count,
                 n13 = fungi_new_york_rg_genus_count,
                 n123 = fungi_all3_genus_count,
                 category = c("new_york_muscat", "c3309", "riparia_gloire"),
                 scaled=TRUE,
                 fill = c( "#009E73","#E69F00","#56B4E9"))

library(extrafont)
require(cowplot)
pdf("rootstock_venn_diagram.pdf", width=6.5, height=3,family="Arial")
plot_grid(grobTree(bacteria_rootstock_venn),grobTree(fungi_rootstock_venn),labels="AUTO")
dev.off()

#Make venn diagrams to show overlapping genera for root_depths and root depths (based on bacteria and fungi). 

library("VennDiagram")

source("root_depth_project_functions.R")

###BACTERIA ROOT DEPTH####

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs
bacteria_ASVs <- read.table("bacteria/dada2_output_exported_grape/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

#Keep only meta data for grape roots
bacteria_meta <- bacteria_meta[colnames(bacteria_ASVs),]

table(colnames(bacteria_ASVs)== rownames(bacteria_meta))
#they match

#Get taxonoomy information
bacteria_taxa <- read.table("bacteria/taxa/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

# Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
bacteria_depth_0_15_samples <- rownames(bacteria_meta)[which(bacteria_meta$root_depth == "depth_0_15")]
bacteria_depth_15_30_samples <- rownames(bacteria_meta)[which(bacteria_meta$root_depth == "depth_15_30")]
bacteria_depth_30_50_samples <- rownames(bacteria_meta)[which(bacteria_meta$root_depth == "depth_30_50")]

depth_0_15_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_depth_0_15_samples]) > 0)]
depth_0_15_genera <- unique(bacteria_taxa_breakdown[depth_0_15_ASVs, "genus"])

depth_15_30_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_depth_15_30_samples]) > 0)]
depth_15_30_genera <- unique(bacteria_taxa_breakdown[depth_15_30_ASVs, "genus"])

depth_30_50_ASVs <- rownames(bacteria_ASVs)[which(rowSums(bacteria_ASVs[, bacteria_depth_30_50_samples]) > 0)]
depth_30_50_genera <- unique(bacteria_taxa_breakdown[depth_30_50_ASVs, "genus"])

bacteria_depth_0_15_genus_count <- length(depth_0_15_genera)
bacteria_depth_15_30_genus_count <- length(depth_15_30_genera)
bacteria_depth_30_50_genus_count <- length(depth_30_50_genera)
bacteria_depth_0_15_depth_15_30_genus_count <- length(which(depth_0_15_genera %in% depth_15_30_genera))
bacteria_depth_15_30_depth_30_50_genus_count <- length(which(depth_15_30_genera %in% depth_30_50_genera))
bacteria_depth_0_15_depth_30_50_genus_count <- length(which(depth_0_15_genera %in% depth_30_50_genera))
bacteria_all3_genus_count_TMP <- depth_0_15_genera[which(depth_0_15_genera %in% depth_30_50_genera)]
bacteria_all3_genus_count <- length(bacteria_all3_genus_count_TMP[which(bacteria_all3_genus_count_TMP %in% depth_15_30_genera)])

# Saved as 7x7 inches
bacteria_root_depth_venn <- draw.triple.venn(area1=bacteria_depth_0_15_genus_count,
                                             area2=bacteria_depth_15_30_genus_count,
                                             area3=bacteria_depth_30_50_genus_count,
                                             n12 = bacteria_depth_0_15_depth_15_30_genus_count,
                                             n23 = bacteria_depth_15_30_depth_30_50_genus_count,
                                             n13 = bacteria_depth_0_15_depth_30_50_genus_count,
                                             n123 = bacteria_all3_genus_count,
                                             category = c("depth_0_15", "depth_15_30", "depth_30_50"),
                                             scaled=TRUE,
                                             fill = c( "#009E73","#E69F00","#56B4E9"))

###FUNGI ROOT DEPTH####

fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs
fungi_ASVs <- read.table("fungi/dada2_output_exported_grape/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

#Keep only meta data for grape roots
fungi_meta <- fungi_meta[colnames(fungi_ASVs),]

table(colnames(fungi_ASVs)== rownames(fungi_meta))
#they match

#Get taxonoomy information
fungi_taxa <- read.table("fungi/taxa/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- qiime2_taxa_breakdown(fungi_taxa)

# Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
fungi_depth_0_15_samples <- rownames(fungi_meta)[which(fungi_meta$root_depth == "depth_0_15")]
fungi_depth_15_30_samples <- rownames(fungi_meta)[which(fungi_meta$root_depth == "depth_15_30")]
fungi_depth_30_50_samples <- rownames(fungi_meta)[which(fungi_meta$root_depth == "depth_30_50")]

depth_0_15_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_depth_0_15_samples]) > 0)]
depth_0_15_genera <- unique(fungi_taxa_breakdown[depth_0_15_ASVs, "genus"])

depth_15_30_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_depth_15_30_samples]) > 0)]
depth_15_30_genera <- unique(fungi_taxa_breakdown[depth_15_30_ASVs, "genus"])

depth_30_50_ASVs <- rownames(fungi_ASVs)[which(rowSums(fungi_ASVs[, fungi_depth_30_50_samples]) > 0)]
depth_30_50_genera <- unique(fungi_taxa_breakdown[depth_30_50_ASVs, "genus"])

fungi_depth_0_15_genus_count <- length(depth_0_15_genera)
fungi_depth_15_30_genus_count <- length(depth_15_30_genera)
fungi_depth_30_50_genus_count <- length(depth_30_50_genera)
fungi_depth_0_15_depth_15_30_genus_count <- length(which(depth_0_15_genera %in% depth_15_30_genera))
fungi_depth_15_30_depth_30_50_genus_count <- length(which(depth_15_30_genera %in% depth_30_50_genera))
fungi_depth_0_15_depth_30_50_genus_count <- length(which(depth_0_15_genera %in% depth_30_50_genera))
fungi_all3_genus_count_TMP <- depth_0_15_genera[which(depth_0_15_genera %in% depth_30_50_genera)]
fungi_all3_genus_count <- length(fungi_all3_genus_count_TMP[which(fungi_all3_genus_count_TMP %in% depth_15_30_genera)])

# Saved as 7x7 inches
fungi_root_depth_venn <- draw.triple.venn(area1=fungi_depth_0_15_genus_count,
                                          area2=fungi_depth_15_30_genus_count,
                                          area3=fungi_depth_30_50_genus_count,
                                          n12 = fungi_depth_0_15_depth_15_30_genus_count,
                                          n23 = fungi_depth_15_30_depth_30_50_genus_count,
                                          n13 = fungi_depth_0_15_depth_30_50_genus_count,
                                          n123 = fungi_all3_genus_count,
                                          category = c("depth_0_15", "depth_15_30", "depth_30_50"),
                                          scaled=TRUE,
                                          fill = c( "#009E73","#E69F00","#56B4E9"))

library(extrafont)
require(cowplot)
pdf("root_depth_venn_diagram.pdf", width=6.5, height=3,family="Arial")
plot_grid(grobTree(bacteria_root_depth_venn),grobTree(fungi_root_depth_venn),labels="AUTO")
dev.off()



