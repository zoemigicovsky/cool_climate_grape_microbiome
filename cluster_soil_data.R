### Cluster soil data to remove redundant features.

rm(list=ls(all=TRUE))

library("factoextra")
library("pheatmap")

# Read in soil data, scale, and get euclidean dist.
soil_data <- read.table("/home/gavin/github_repos/root_depth/soil_dat.tsv",
                        header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

soil_data_scaled <- scale(soil_data, center = TRUE, scale = TRUE)

soil_data_scaled_dist <- dist(soil_data_scaled, method = "euclidean")

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
                                    cutree_cols=6)

# Identified features to retain: c("fe", "al", "mg_sat", "na", "cu", "ca")
