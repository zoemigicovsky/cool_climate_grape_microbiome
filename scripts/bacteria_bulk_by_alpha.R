rm(list=ls(all=TRUE))

library(cowplot)
library(ggplot2)
library(ggbeeswarm)
library(Hmisc)

setwd("/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/")
source("/home/gavin/github_repos/root_depth/root_depth_project_functions.R")

bacteria_richness <- read.table("bacteria/diversity_no_chl_no_mit/observed_otus_vector_exported/alpha-diversity.tsv", header=TRUE, sep="\t", row.names=1)
bacteria_evenness <- read.table("bacteria/diversity_no_chl_no_mit/evenness_vector_exported/alpha-diversity.tsv", header=TRUE, sep="\t", row.names=1)
bacteria_phylo_d <- read.table("bacteria/diversity_no_chl_no_mit/faith_pd_vector_exported/alpha-diversity.tsv", header=TRUE, sep="\t", row.names=1)

fungi_richness <- read.table("fungi/diversity_no_chl_no_mit/observed_otus_vector_exported/alpha-diversity.tsv", header=TRUE, sep="\t", row.names=1)
fungi_evenness <- read.table("fungi/diversity_no_chl_no_mit/evenness_vector_exported/alpha-diversity.tsv", header=TRUE, sep="\t", row.names=1)
fungi_phylo_d <- read.table("fungi/diversity_no_chl_no_mit/faith_pd_vector_exported/alpha-diversity.tsv", header=TRUE, sep="\t", row.names=1)

fungi_abun <- read.table("fungi/dada2_output_exported/feature-table_w_tax.txt", header=TRUE, row.names=1, comment.char="", skip=1, sep="\t")
fungi_abun <- fungi_abun[, -which(colnames(fungi_abun) == "taxonomy")]

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

bacteria_meta_soil <- bacteria_meta[which(bacteria_meta$tissue == "soil"), c("rootstock", "root_depth")]
fungi_meta_soil <- fungi_meta[which(fungi_meta$tissue == "soil"), c("rootstock", "root_depth")]

bacteria_meta_soil <- bacteria_meta_soil[which(rownames(bacteria_meta_soil) %in% rownames(bacteria_richness)), ]
fungi_meta_soil <- fungi_meta_soil[which(rownames(fungi_meta_soil) %in% rownames(fungi_richness)), ]

bacteria_alpha <- bacteria_meta_soil
bacteria_alpha$richness <- bacteria_richness[rownames(bacteria_alpha), 1]
bacteria_alpha$evenness <- bacteria_evenness[rownames(bacteria_alpha), 1]
bacteria_alpha$phylo_d <- bacteria_phylo_d[rownames(bacteria_alpha), 1]

bacteria_root_depth_richness_lm <- lm(formula = richness ~ root_depth, data = bacteria_alpha)
bacteria_root_depth_evenness_lm <- lm(formula = evenness ~ root_depth, data = bacteria_alpha)
bacteria_root_depth_phylo_d_lm <- lm(formula = phylo_d ~ root_depth, data = bacteria_alpha)

summary(bacteria_root_depth_richness_lm)
summary(bacteria_root_depth_evenness_lm)
summary(bacteria_root_depth_phylo_d_lm)

fungi_alpha <- fungi_meta_soil
fungi_alpha$richness <- fungi_richness[rownames(fungi_alpha), 1]
fungi_alpha$evenness <- fungi_evenness[rownames(fungi_alpha), 1]
fungi_alpha$phylo_d <- fungi_phylo_d[rownames(fungi_alpha), 1]

fungi_root_depth_richness_lm <- lm(formula = richness ~ root_depth, data = fungi_alpha)
fungi_root_depth_evenness_lm <- lm(formula = evenness ~ root_depth, data = fungi_alpha)
fungi_root_depth_phylo_d_lm <- lm(formula = phylo_d ~ root_depth, data = fungi_alpha)

summary(fungi_root_depth_richness_lm)
summary(fungi_root_depth_evenness_lm)
summary(fungi_root_depth_phylo_d_lm)


# Clean up metadata.
colnames(bacteria_alpha)[which(colnames(bacteria_alpha) == "rootstock")] <- "Rootstock"
bacteria_alpha[which(bacteria_alpha$Rootstock == "c3309"), "Rootstock"] <- "C3309"
bacteria_alpha[which(bacteria_alpha$Rootstock == "new_york_muscat"), "Rootstock"] <- "Ungrafted"
bacteria_alpha[which(bacteria_alpha$Rootstock == "riparia_gloire"), "Rootstock"] <- "Riparia Gloire"

bacteria_alpha$Rootstock <- factor(bacteria_alpha$Rootstock, levels=c("Ungrafted", "C3309", "Riparia Gloire"))

bacteria_alpha[which(bacteria_alpha$root_depth == "depth_0_15"), "root_depth"] <- "0-15 cm"
bacteria_alpha[which(bacteria_alpha$root_depth == "depth_15_30"), "root_depth"] <- "15-30 cm"
bacteria_alpha[which(bacteria_alpha$root_depth == "depth_30_50"), "root_depth"] <- "30-50 cm"

colnames(fungi_alpha)[which(colnames(fungi_alpha) == "rootstock")] <- "Rootstock"
fungi_alpha[which(fungi_alpha$Rootstock == "c3309"), "Rootstock"] <- "C3309"
fungi_alpha[which(fungi_alpha$Rootstock == "new_york_muscat"), "Rootstock"] <- "Ungrafted"
fungi_alpha[which(fungi_alpha$Rootstock == "riparia_gloire"), "Rootstock"] <- "Riparia Gloire"

fungi_alpha$Rootstock <- factor(fungi_alpha$Rootstock, levels=c("Ungrafted", "C3309", "Riparia Gloire"))

fungi_alpha[which(fungi_alpha$root_depth == "depth_0_15"), "root_depth"] <- "0-15 cm"
fungi_alpha[which(fungi_alpha$root_depth == "depth_15_30"), "root_depth"] <- "15-30 cm"
fungi_alpha[which(fungi_alpha$root_depth == "depth_30_50"), "root_depth"] <- "30-50 cm"

# Make plots.
legend_text_size <- 8

bacteria_richness_rsq <- round(summary(bacteria_root_depth_richness_lm)$adj.r.squared, digits=3)
bacteria_richness_f <- summary(bacteria_root_depth_richness_lm)$fstatistic
bacteria_richness_p <- round(as.numeric(pf(bacteria_richness_f[1], bacteria_richness_f[2], bacteria_richness_f[3], lower.tail=F)), digits=3)

#bacteria_richness_plot_title <- bquote("Bacteria richness (Adj." ~ R^2~"=" ~ .(bacteria_richness_rsq) ~ ", P=" ~ .(bacteria_richness_p) ~ ")")
bacteria_richness_plot_title <- bquote("Bacteria richness")

bacteria_richness <- ggplot(aes(colour=Rootstock, x=root_depth, y=richness), data=bacteria_alpha) +
                                geom_quasirandom(size=5, width=0.1) +
                                stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                                             geom="pointrange", color="black",  alpha=0.4) +
                              ylab("Richness") +
                              xlab("Root depth") +
                              ggtitle(bacteria_richness_plot_title) +
  theme(legend.position=c(0.02, 0.17),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size))



bacteria_evenness_rsq <- round(summary(bacteria_root_depth_evenness_lm)$adj.r.squared, digits=3)
bacteria_evenness_f <- summary(bacteria_root_depth_evenness_lm)$fstatistic
bacteria_evenness_p <- round(as.numeric(pf(bacteria_evenness_f[1], bacteria_evenness_f[2], bacteria_evenness_f[3], lower.tail=F)), digits=3)

#bacteria_evenness_plot_title <- bquote("Bacteria evenness (Adj." ~ R^2~"=" ~ .(bacteria_evenness_rsq) ~ ", P=" ~ .(bacteria_evenness_p) ~ ")")
bacteria_evenness_plot_title <- bquote("Bacteria evenness")

bacteria_evenness <- ggplot(aes(colour=Rootstock, x=root_depth, y=evenness), data=bacteria_alpha) +
  geom_quasirandom(size=5, width=0.1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="black",  alpha=0.4) +
  ylab("Evenness") +
  xlab("Root depth") +
  ggtitle(bacteria_evenness_plot_title) +
  theme(legend.position=c(0.65, 0.17),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size))


bacteria_phylo_d_rsq <- round(summary(bacteria_root_depth_phylo_d_lm)$adj.r.squared, digits=3)
bacteria_phylo_d_f <- summary(bacteria_root_depth_phylo_d_lm)$fstatistic
bacteria_phylo_d_p <- round(as.numeric(pf(bacteria_phylo_d_f[1], bacteria_phylo_d_f[2], bacteria_phylo_d_f[3], lower.tail=F)), digits=3)

#bacteria_phylo_d_plot_title <- bquote("Bacteria PD (Adj." ~ R^2~"=" ~ .(bacteria_phylo_d_rsq) ~ ", P=" ~ .(bacteria_phylo_d_p) ~ ")")
bacteria_phylo_d_plot_title <- bquote("Bacteria PD")

bacteria_phylo_d <- ggplot(aes(colour=Rootstock, x=root_depth, y=phylo_d), data=bacteria_alpha) +
  geom_quasirandom(size=5, width=0.1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="black",  alpha=0.4) +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("Root depth") +
  ggtitle(bacteria_phylo_d_plot_title) +
  theme(legend.position=c(0.37, 0.17),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size))

fungi_richness_rsq <- round(summary(fungi_root_depth_richness_lm)$adj.r.squared, digits=3)
fungi_richness_f <- summary(fungi_root_depth_richness_lm)$fstatistic
fungi_richness_p <- round(as.numeric(pf(fungi_richness_f[1], fungi_richness_f[2], fungi_richness_f[3], lower.tail=F)), digits=3)

#fungi_richness_plot_title <- bquote("Fungi richness (Adj." ~ R^2~"=" ~ .(fungi_richness_rsq) ~ ", P=" ~ .(fungi_richness_p) ~ ")")
fungi_richness_plot_title <- bquote("Fungi richness")

fungi_richness <- ggplot(aes(colour=Rootstock, x=root_depth, y=richness), data=fungi_alpha) +
  geom_quasirandom(size=5, width=0.1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="black",  alpha=0.4) +
  ylab("Richness") +
  xlab("Root depth") +
  ggtitle(fungi_richness_plot_title) +
  theme(legend.position=c(0.02, 0.17),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size))



fungi_evenness_rsq <- round(summary(fungi_root_depth_evenness_lm)$adj.r.squared, digits=3)
fungi_evenness_f <- summary(fungi_root_depth_evenness_lm)$fstatistic
fungi_evenness_p <- round(as.numeric(pf(fungi_evenness_f[1], fungi_evenness_f[2], fungi_evenness_f[3], lower.tail=F)), digits=3)

#fungi_evenness_plot_title <- bquote("Fungi evenness (Adj." ~ R^2~"=" ~ .(fungi_evenness_rsq) ~ ", P=" ~ .(fungi_evenness_p) ~ ")")
fungi_evenness_plot_title <- bquote("Fungi evenness")

fungi_evenness <- ggplot(aes(colour=Rootstock, x=root_depth, y=evenness), data=fungi_alpha) +
  geom_quasirandom(size=5, width=0.1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="black",  alpha=0.4) +
  ylab("Evenness") +
  xlab("Root depth") +
  ggtitle(fungi_evenness_plot_title) +
  theme(legend.position=c(0.65, 0.17),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size))


fungi_phylo_d_rsq <- round(summary(fungi_root_depth_phylo_d_lm)$adj.r.squared, digits=3)
fungi_phylo_d_f <- summary(fungi_root_depth_phylo_d_lm)$fstatistic
fungi_phylo_d_p <- round(as.numeric(pf(fungi_phylo_d_f[1], fungi_phylo_d_f[2], fungi_phylo_d_f[3], lower.tail=F)), digits=3)

#fungi_phylo_d_plot_title <- bquote("Fungi PD (Adj." ~ R^2~"=" ~ .(fungi_phylo_d_rsq) ~ ", P=" ~ .(fungi_phylo_d_p) ~ ")")
fungi_phylo_d_plot_title <- bquote("Fungi PD")

fungi_phylo_d <- ggplot(aes(colour=Rootstock, x=root_depth, y=phylo_d), data=fungi_alpha) +
  geom_quasirandom(size=5, width=0.1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="black",  alpha=0.4) +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("Root depth") +
  ggtitle(fungi_phylo_d_plot_title) +
  theme(legend.position=c(0.02, 0.17),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size))

plot_grid(bacteria_richness, bacteria_evenness, bacteria_phylo_d,
          fungi_richness, fungi_evenness, fungi_phylo_d,
          labels=c("A", "B", "C", "D", "E", "F"),
          nrow=2,
          ncol=3)
