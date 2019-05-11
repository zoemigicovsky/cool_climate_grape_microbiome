rm(list=ls(all=TRUE))

library("ggvegan")
library("vegan")
library("reshape2")
library("cowplot")

setwd("/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/")

# Read in soil data.
soil_data <- read.table("/home/gavin/github_repos/root_depth/soil_dat.tsv",
                        header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

### First process 16S data.
# Read in BIOM.
biom_16S <- read.table("bacteria/dada2_output_exported/feature-table_w_tax.txt", header=T, sep="\t", row.names=1, skip=1, comment.char="")
biom_16S <- biom_16S[ , -which(colnames(biom_16S) == "taxonomy")]

# Read in sample mapfile.
map_16S <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                      header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

# Subset tables to overlapping samples across soil data and BIOM.
overlapping_samples_16S <- rownames(soil_data)[which(rownames(soil_data) %in% colnames(biom_16S))]

biom_16S <- biom_16S[, overlapping_samples_16S]
biom_16S <- biom_16S[-which(rowSums(biom_16S) == 0), ]
soil_16S <- soil_data[overlapping_samples_16S, ]
map_16S <- map_16S[overlapping_samples_16S, ]

### Transpose OTU table and convert to rel. abun.
biom_16S_t <- data.frame(t(biom_16S), check.names=FALSE)
biom_16S_t_relabun <- data.frame(sweep(biom_16S_t, 1, rowSums(biom_16S_t), '/'), check.names=FALSE) * 100

cca_16S_full <- cca(biom_16S_t_relabun ~ ., soil_16S)

cca_16S_fortify <- ggvegan:::fortify.cca(cca_16S_full)

cca_16S_fortify_sites <- cca_16S_fortify[which(cca_16S_fortify$Score == "sites"), ]

cca_16S_fortify_sites$rootstock <- NA
cca_16S_fortify_sites[, "rootstock"] <- map_16S[as.character(cca_16S_fortify_sites$Label), "rootstock"]
cca_16S_fortify_sites$rootstock <- as.factor(cca_16S_fortify_sites$rootstock)

cca_16S_fortify_sites$root_depth <- NA
cca_16S_fortify_sites[, "root_depth"] <- map_16S[as.character(cca_16S_fortify_sites$Label), "root_depth"]
cca_16S_fortify_sites$root_depth <- as.factor(cca_16S_fortify_sites$root_depth)

# Get full CCA plot
cca_16S_plot <- ggvegan:::autoplot.cca(cca_16S_fortify, layers=c("biplot")) +
                                       geom_point(data=cca_16S_fortify_sites, aes(x=CCA1, y=CCA2, col=rootstock, shape=root_depth), size=2) +
                                       scale_colour_manual(name = "Rootstock", values=c("orange", "blue", "dark green")) +
                                       scale_shape_discrete(name = "Depth") +
                                       ggtitle("Bacteria (16S)") +
                                       theme_bw() +
                                       theme(plot.title = element_text(face="bold")) #+
                                       #guides(colour=FALSE, shape=FALSE)

### Do same for ITS data.
# Read in BIOM.
biom_ITS <- read.table("fungi/dada2_output_exported/feature-table_w_tax.txt", header=T, sep="\t", row.names=1, skip=1, comment.char="")
biom_ITS <- biom_ITS[ , -which(colnames(biom_ITS) == "taxonomy")]

# Read in sample mapfile.
map_ITS <- read.table("fungi/root_depth_fungi_metadata.tsv",
                      header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

# Subset tables to overlapping samples across soil data and BIOM.
overlapping_samples_ITS <- rownames(soil_data)[which(rownames(soil_data) %in% colnames(biom_ITS))]

biom_ITS <- biom_ITS[, overlapping_samples_ITS]
biom_ITS <- biom_ITS[-which(rowSums(biom_ITS) == 0), ]
soil_ITS <- soil_data[overlapping_samples_ITS, ]
map_ITS <- map_ITS[overlapping_samples_ITS, ]

### Transpose OTU table and convert to rel. abun.
biom_ITS_t <- data.frame(t(biom_ITS), check.names=FALSE)
biom_ITS_t_relabun <- data.frame(sweep(biom_ITS_t, 1, rowSums(biom_ITS_t), '/'), check.names=FALSE) * 100

### Need to cluster soil data based on Spearman correlations and collapse to a small number of variables.
cca_ITS_full <- cca(biom_ITS_t_relabun ~ ., soil_ITS)

cca_ITS_fortify <- ggvegan:::fortify.cca(cca_ITS_full)

cca_ITS_fortify_sites <- cca_ITS_fortify[which(cca_ITS_fortify$Score == "sites"), ]

cca_ITS_fortify_sites$rootstock <- NA
cca_ITS_fortify_sites[, "rootstock"] <- map_ITS[as.character(cca_ITS_fortify_sites$Label), "rootstock"]
cca_ITS_fortify_sites$rootstock <- as.factor(cca_ITS_fortify_sites$rootstock)

cca_ITS_fortify_sites$root_depth <- NA
cca_ITS_fortify_sites[, "root_depth"] <- map_ITS[as.character(cca_ITS_fortify_sites$Label), "root_depth"]
cca_ITS_fortify_sites$root_depth <- as.factor(cca_ITS_fortify_sites$root_depth)

# Get full CCA plot.
cca_ITS_plot <- ggvegan:::autoplot.cca(cca_ITS_fortify, layers=c("biplot")) +
  geom_point(data=cca_ITS_fortify_sites, aes(x=CCA1, y=CCA2, col=rootstock, shape=root_depth), size=2) +
  scale_colour_manual(name = "Rootstock", values=c("orange", "blue", "dark green")) +
  scale_shape_discrete(name = "Depth") +
  ggtitle("Fungi (ITS)") +
  theme_bw() +
  theme(plot.title = element_text(face="bold")) #+
  #guides(colour=FALSE, shape=FALSE)

### Plot figure - I couldn't get this looking pretty so left it for now.
# Retain only legend for downstream plotting.
#   cca_16S_plot_TMP <- ggvegan:::autoplot.cca(cca_16S_fortify, layers=c("biplot")) +
#   geom_point(data=cca_16S_fortify_sites, aes(x=CCA1, y=CCA2, col=rootstock, shape=root_depth), size=2) +
#   scale_colour_manual(name = "Rootstock", values=c("orange", "blue", "dark green")) +
#   scale_shape_discrete(name = "Depth")
# 
# grobs <- ggplotGrob(cca_16S_plot_TMP)$grobs
# cca_plot_legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
#   plot_grid(cca_16S_plot, cca_ITS_plot, cca_plot_legend, labels = c('A', 'B', ''), rel_widths = c(1, 1, 1), nrow = 1)
