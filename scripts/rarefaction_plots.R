### Commands to make rarefaction plots for both bacteria and fungi samples.

rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(randomcoloR)
library(reshape2)
library(stringr)

# Change this line for different systems:
path2repo <- "/home/gavin/github_repos/root_depth"

setwd(path2repo)

bacteria_rarefaction <- read.table("data/rarefaction_data/bacteria_observed_otus.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE)

bacteria_rarefaction_melt <- melt(bacteria_rarefaction)
bacteria_rarefaction_melt$variable <- gsub("_iter.*$", "", bacteria_rarefaction_melt$variable)

bacteria_rarefaction_mean <- aggregate(formula = value ~ variable + sample.id,
                                       data = bacteria_rarefaction_melt,
                                       FUN = mean)

bacteria_rarefaction_mean$depth <- as.numeric(gsub("depth.", "", bacteria_rarefaction_mean$variable))

bacteria_n <- length(bacteria_rarefaction_mean$sample.id[-which(duplicated(bacteria_rarefaction_mean$sample.id))])
bacteria_palette <- distinctColorPalette(bacteria_n)

bacteria_rarefaction_plot <- ggplot(data = bacteria_rarefaction_mean, aes(x = depth, y = value, group = sample.id, colour = sample.id)) +
                                    geom_line() +
                                    geom_point() +
                                    geom_vline(xintercept = 4000, linetype = "dashed", color = "dark green") +
                                    theme_bw() +
                                    xlab("Read depth") +
                                    ylab("Number of ASVs") +
                                    scale_colour_manual(values = bacteria_palette) +
                                    theme(legend.position = 'none') +
                                    ggtitle("Bacteria")



fungi_rarefaction <- read.table("data/rarefaction_data/fungi_observed_otus.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE)

fungi_rarefaction_melt <- melt(fungi_rarefaction)
fungi_rarefaction_melt$variable <- gsub("_iter.*$", "", fungi_rarefaction_melt$variable)

fungi_rarefaction_mean <- aggregate(formula = value ~ variable + sample.id,
                                       data = fungi_rarefaction_melt,
                                       FUN = mean)

fungi_rarefaction_mean$depth <- as.numeric(gsub("depth.", "", fungi_rarefaction_mean$variable))

fungi_n <- length(fungi_rarefaction_mean$sample.id[-which(duplicated(fungi_rarefaction_mean$sample.id))])
fungi_palette <- distinctColorPalette(fungi_n)

fungi_rarefaction_plot <- ggplot(data = fungi_rarefaction_mean, aes(x = depth, y = value, group = sample.id, colour = sample.id)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 4000, linetype = "dashed", color = "dark green") +
  theme_bw() +
  xlab("Read depth") +
  ylab("Number of ASVs") +
  scale_colour_manual(values = fungi_palette) +
  theme(legend.position = 'none') +
  ggtitle("Fungi")

rarefaction_plot <- plot_grid(bacteria_rarefaction_plot,
                              fungi_rarefaction_plot,
                              labels = c('a', 'b'))

ggsave(filename = "figures/rarefaction_plots.pdf",
       plot = rarefaction_plot,
       device = "pdf", dpi = 600, width = 8, heigh = 4, units = "in")

ggsave(filename = "figures/rarefaction_plots.png",
       plot = rarefaction_plot,
       device = "png", dpi = 300, width = 8, heigh = 4, units = "in")

