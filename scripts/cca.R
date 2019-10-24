rm(list=ls(all=TRUE))

library("ggvegan")
library("vegan")
library("reshape2")
library("cowplot")

setwd("/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/")

# Function to return CCA plot.
cca_rootstock_rootdepth <- function(biom_tab, soil_info, metadata, plot_title, no_legend=TRUE) {
  
  cca_full <- cca(biom_tab ~ ., soil_info)
  
  cca_fortify <- ggvegan:::fortify.cca(cca_full)
  
  cca_fortify_sites <- cca_fortify[which(cca_fortify$Score == "sites"), ]
  
  cca_fortify_sites$rootstock <- NA
  cca_fortify_sites[, "rootstock"] <- metadata[as.character(cca_fortify_sites$Label), "rootstock"]
  cca_fortify_sites$rootstock <- as.factor(cca_fortify_sites$rootstock)
  
  cca_fortify_sites$root_depth <- NA
  cca_fortify_sites[, "root_depth"] <- metadata[as.character(cca_fortify_sites$Label), "root_depth"]
  cca_fortify_sites$root_depth <- as.factor(cca_fortify_sites$root_depth)
  
  # Get full CCA plot
  cca_plot <- ggvegan:::autoplot.cca(cca_fortify, layers=c("biplot")) +
    geom_point(data=cca_fortify_sites, aes(x=CCA1, y=CCA2, col=rootstock, shape=root_depth), size=2) +
    scale_colour_manual(name = "Rootstock", values=c("orange", "blue", "dark green")) +
    scale_shape_discrete(name = "Depth") +
    ggtitle(plot_title) +
    theme_bw() +
    theme(plot.title = element_text(face="bold"))
  
  if(no_legend) {
    cca_plot <- cca_plot + theme(legend.position = "none") 
  }
  
  return(cca_plot)
}

# Read in soil data.
soil_data <- read.table("/home/gavin/github_repos/root_depth/soil_dat.tsv",
                        header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)

# Clean up soil data names.
soil_naming <- read.table("/home/gavin/github_repos/root_depth/soil_factor_names.tsv",
                          header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", row.names=1)
colnames(soil_data) <- soil_naming[colnames(soil_data), ]

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

# Get CCA plots.
soil_factors <- c("Fe", "Al", "Mg_sat", "Na", "Mg", "OM")
soil_factors_w_K_sat <- c("Fe", "Al", "Mg_sat", "Na", "Mg", "OM", "K_sat")
soil_factors_w_P2O5 <- c("Fe", "Al", "Mg_sat", "Na", "Mg", "OM", "P2O5")

cca_16S_plot <- cca_rootstock_rootdepth(biom_tab=biom_16S_t_relabun,
                                            soil_info=soil_16S[, soil_factors],
                                            metadata=map_16S,
                                            plot_title="Bacteria (16S)")

cca_16S_plot_w_K_sat <- cca_rootstock_rootdepth(biom_tab=biom_16S_t_relabun,
                                        soil_info=soil_16S[, soil_factors_w_K_sat],
                                        metadata=map_16S,
                                        plot_title="Bacteria (16S) - with Ksat")

cca_16S_plot_w_P2O5 <- cca_rootstock_rootdepth(biom_tab=biom_16S_t_relabun,
                                               soil_info=soil_16S[, soil_factors_w_P2O5],
                                               metadata=map_16S,
                                               plot_title="Bacteria (16S) - with P2O5")


cca_ITS_plot <- cca_rootstock_rootdepth(biom_tab=biom_ITS_t_relabun,
                                        soil_info=soil_ITS[, soil_factors],
                                        metadata=map_ITS,
                                        plot_title="Fungi (ITS)")

cca_ITS_plot_w_K_sat <- cca_rootstock_rootdepth(biom_tab=biom_ITS_t_relabun,
                                                soil_info=soil_ITS[, soil_factors_w_K_sat],
                                                metadata=map_ITS,
                                                plot_title="Fungi (ITS) - with Ksat")

cca_ITS_plot_w_P2O5 <- cca_rootstock_rootdepth(biom_tab=biom_ITS_t_relabun,
                                               soil_info=soil_ITS[, soil_factors_w_P2O5],
                                               metadata=map_ITS,
                                               plot_title="Fungi (ITS) - with P2O5")


# Parse legend for plotting.
cca_16S_plot_TMP <- cca_rootstock_rootdepth(biom_tab=biom_16S_t_relabun,
                                            soil_info=soil_16S[, soil_factors],
                                            metadata=map_16S,
                                            plot_title="Bacteria (16S)",
                                            no_legend=FALSE)

grobs <- ggplotGrob(cca_16S_plot_TMP)$grobs
cca_plot_legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

plot_grid(cca_16S_plot, cca_16S_plot_w_K_sat, cca_16S_plot_w_P2O5, cca_plot_legend, ncol=2, labels=c("", "", "", ""))

plot_grid(cca_ITS_plot, cca_ITS_plot_w_K_sat, cca_ITS_plot_w_P2O5, cca_plot_legend, ncol=2, labels=c("", "", "", ""))
