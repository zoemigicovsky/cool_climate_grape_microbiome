### Cluster soil data to remove redundant features.

rm(list=ls(all=TRUE))

library("factoextra")
library("pheatmap")

# Read in soil data, scale, and get euclidean dist.
soil_data <- read.table("/home/gavin/github_repos/root_depth/soil_dat.tsv",
                        header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

soil_data_scaled <- scale(soil_data, center = TRUE, scale = TRUE)

soil_data_scaled_dist <- dist(soil_data_scaled, method = "euclidean")


# Make metadata table with rootdepth and rootstock.
map_16S <- read.table("/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/bacteria/root_depth_bacteria_metadata.tsv",
                      header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

map_16S_subset <- map_16S[rownames(soil_data), c("rootstock", "root_depth")]

# Make elbow plot to help choose number of clusters to retain.
fviz_nbclust(x=soil_data_scaled,
             diss=soil_data_scaled_dist,
             hcut,
             method = "wss",
             k.max = 20)

# Cluster soil data based on euclidean distances of scaled data.
soil_data_scaled_pheatmap <- pheatmap(soil_data_scaled,
                                    clustering_distance_cols = "euclidean",
                                    clustering_method = "complete",
                                    cutree_cols=6,
                                    annotation_row=map_16S_subset)

# Identified features to retain: c("fe", "al", "mg_sat", "na", "cu", "ca")
