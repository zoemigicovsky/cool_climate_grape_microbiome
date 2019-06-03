#Figure 1 will be a comparison across grape roots/soil/cover crops

#Panel A will be a venn diagram (bacteria)

#Panel B will be stacked bar charts showing main differences (bacteria)

#Panel C will be a venn diagram (fungi)

#Panel D will be stacked bar charts showing main differences (fungi)

###Venn diagram - Bacteria####

# Read in files and determine genera overlapping in cover, soil, and root for bacteria datasets.

rm(list=ls(all=TRUE))

library(extrafont)
library(cowplot)
library(vegan)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(inlmisc)

# Change these two lines for different systems:
path2repo <- "/home/gavin/github_repos/root_depth"
working_dir <- "/home/gavin/gavin_backup/projects/zoe_microbiome/data/root_depth/"

source(paste(path2repo, "root_depth_project_functions.R", sep="/"))
setwd(working_dir)

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

bacteria_meta$group <- bacteria_meta$tissue
bacteria_meta$group[which(bacteria_meta$species == "cover_crop")] <- "root (cover)"

bacteria_ASVs <- read.table("bacteria/dada2_output_exported/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

# Get min per-sample read depth.
min(colSums(bacteria_ASVs))
# 4564

# Rarefy table. Save this output table as RDS file to avoid re-running this step.
# Output RDS filepath:
bacteria_ASVs_4564subsample_RDS <- paste(path2repo, "intermediate_RDS/bacteria_ASVs_4564subsample.rds", sep="/")

# # Code that should not be re-run every time:
# set.seed(141)
# bacteria_ASVs_4564subsample <- data.frame(t(rrarefy(x = t(bacteria_ASVs), sample=4564)), check.names=FALSE)
# bacteria_ASVs_4564subsample <- bacteria_ASVs_4564subsample[-which(rowSums(bacteria_ASVs_4564subsample) == 0), ]
# saveRDS(object = bacteria_ASVs_4564subsample, file = bacteria_ASVs_4564subsample_RDS)

# Read in RDS:
bacteria_ASVs_4564subsample <- readRDS(bacteria_ASVs_4564subsample_RDS)


bacteria_taxa <- read.table("bacteria/taxa/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

bacteria_sampletype_venn <- threeWayVennPercent(metadata=bacteria_meta,
                                            meta_col="group",
                                            asv_abun=bacteria_ASVs_4564subsample,
                                            taxa_df=bacteria_taxa_breakdown,
                                            meta_cat=c("root", "root (cover)", "soil"),
                                            labels=c("Grape roots", "Cover crop\nroots", "Grape soil"),
                                            colours=c("#009E73", "#E69F00", "#56B4E9"))

#### Run the same commands, but for fungi:

fungi_ASVs <- read.table("fungi/dada2_output_exported/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

# Get min per-sample read depth.
min(colSums(fungi_ASVs))
# 1176

# Rarefy table. Save this output table as RDS file to avoid re-running this step.
# Output RDS filepath:
fungi_ASVs_1176subsample_RDS <- paste(path2repo, "intermediate_RDS/fungi_ASVs_1176subsample.rds", sep="/")

# Code that should not be re-run every time:
# set.seed(151)
# fungi_ASVs_1176subsample <- data.frame(t(rrarefy(x = t(fungi_ASVs), sample=1176)), check.names=FALSE)
# fungi_ASVs_1176subsample <- fungi_ASVs_1176subsample[-which(rowSums(fungi_ASVs_1176subsample) == 0), ]
# saveRDS(object = fungi_ASVs_1176subsample, file = fungi_ASVs_1176subsample_RDS)

# Read in RDS:
fungi_ASVs_1176subsample <- readRDS(fungi_ASVs_1176subsample_RDS)

fungi_taxa <- read.table("fungi/taxa/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_meta$group <- fungi_meta$tissue
fungi_meta$group[which(fungi_meta$species == "cover_crop")] <- "root (cover)"

# Need to subset metadata to only samples in BIOM table.
fungi_meta <- fungi_meta[which(rownames(fungi_meta) %in% colnames(fungi_ASVs_1176subsample)), ]

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

fungi_sampletype_venn <- threeWayVennPercent(metadata=fungi_meta,
                                                meta_col="group",
                                                asv_abun=fungi_ASVs_1176subsample,
                                                taxa_df=fungi_taxa_breakdown,
                                                meta_cat=c("root", "root (cover)", "soil"),
                                                labels=c("Grape roots", "Cover crop\nroots", "Grape soil"),
                                                colours=c("#009E73", "#E69F00", "#56B4E9"))


###Stacked bar - Bacteria####
bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

bacteria_meta$group <- bacteria_meta$tissue
bacteria_meta$group[which(bacteria_meta$species == "cover_crop")] <- "root (cover)"

bacteria_ASVs <- read.table("bacteria/dada2_output_exported/feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

bacteria_taxa <- read.table("bacteria/taxa/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

#Convert to relative abundance (i.e. sum to 100%) - the "sweep" function is helpful for this.
bacteria_asv_abun <- data.frame(sweep(bacteria_ASVs, 2, colSums(bacteria_ASVs), '/')) * 100

#Next aggregate ASV abundances by genus.

bacteria_asv_abun$genus <- bacteria_taxa_breakdown[rownames(bacteria_asv_abun), "genus"]
bacteria_asv_abun_genus_sum <- aggregate(. ~ genus, data=bacteria_asv_abun, FUN=sum)

#We need to identify rare genera and collapse them into the "Other" category. 

#Reorder based on abundance
bacteria_asv_abun_genus_sum <- bacteria_asv_abun_genus_sum[order(rowSums(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)]),decreasing=T),]
#Keep the top 15 but rename everything else starting at the 16th most abundant as Other 
bacteria_asv_abun_genus_sum[16:nrow(bacteria_asv_abun_genus_sum),1] <- "Other"

#Now melt this table
bacteria_asv_abun_relab_genus_sum_melt <- melt(bacteria_asv_abun_genus_sum)
#Join in meta_data and change sample names 
bacteria_meta <- bacteria_meta %>% mutate(variable=as.character(rownames(bacteria_meta)))
bacteria_asv_abun_relab_genus_sum_melt <-  bacteria_asv_abun_relab_genus_sum_melt %>% inner_join(bacteria_meta)

#reorder samples based on group
bacteria_asv_abun_relab_genus_sum_melt$group <- gsub("root (cover)", "cover", bacteria_asv_abun_relab_genus_sum_melt$group, fixed=TRUE)

bacteria_asv_abun_relab_genus_sum_melt$variable  <- fct_relevel(bacteria_asv_abun_relab_genus_sum_melt$variable ,"e98", "e99")

bacteria_asv_abun_relab_genus_sum_melt$variable  <- fct_relevel(bacteria_asv_abun_relab_genus_sum_melt$variable ,"e104", "e105","e106", "e107","e108", "e109", "e110","e111", "e112", after= Inf)


#get custom colour palette
colour_count = length(unique(bacteria_asv_abun_relab_genus_sum_melt$genus))
cols <- inlmisc::GetColors(n = 16, scheme = "discrete rainbow")
cols[colour_count] <- "darkgrey"
my_palette <- as.vector(cols[1:16])

# Clean up genus labels:
bacteria_asv_abun_relab_genus_sum_melt$genus <- as.character(bacteria_asv_abun_relab_genus_sum_melt$genus)
bacteria_asv_abun_relab_genus_sum_melt$genus <- gsub("D_0__Bacteria;D_1__", "", bacteria_asv_abun_relab_genus_sum_melt$genus)
bacteria_asv_abun_relab_genus_sum_melt$genus <- gsub(";D_.__", "_", bacteria_asv_abun_relab_genus_sum_melt$genus)

#reorder genus based on abudance (low to high, across all samples)
bacteria_asv_abun_relab_genus_sum_melt$genus <- fct_reorder(bacteria_asv_abun_relab_genus_sum_melt$genus, bacteria_asv_abun_relab_genus_sum_melt$value, sum)
bacteria_asv_abun_relab_genus_sum_melt$genus <- fct_relevel(bacteria_asv_abun_relab_genus_sum_melt$genus, "Other", after = Inf)




bacteria_x_col <- rep(NA, length(levels(bacteria_asv_abun_relab_genus_sum_melt$variable)))
names(bacteria_x_col) <- levels(bacteria_asv_abun_relab_genus_sum_melt$variable)
bacteria_x_col[bacteria_meta[which(bacteria_meta$group == "root"), "variable"]] <- "#009E73"
bacteria_x_col[bacteria_meta[which(bacteria_meta$group == "root (cover)"), "variable"]] <- "#E69F00"
bacteria_x_col[bacteria_meta[which(bacteria_meta$group == "soil"), "variable"]] <- "#56B4E9"

#plot
bacteria_stacked <- ggplot(bacteria_asv_abun_relab_genus_sum_melt, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("Bacteria")+
  theme(legend.position="bottom",
        legend.key.size = unit(1, "cm"),
        legend.title=element_blank(),
        axis.text.x = element_text(size=20,
                                   angle = 45,
                                   hjust = 1,
                                   colour=bacteria_x_col),
        axis.text.y = element_text(size=20),
        legend.text=element_text(size=17),
        axis.title=element_text(size=20),
        title=element_text(size=20)) +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  scale_fill_manual(values = my_palette)

###Stacked bar - Fungi####
fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_meta$group <- fungi_meta$tissue
fungi_meta$group[which(fungi_meta$species == "cover_crop")] <- "root (cover)"

fungi_ASVs <- read.table("fungi/dada2_output_exported/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

fungi_taxa <- read.table("fungi/taxa/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

#Convert to relative abundance (i.e. sum to 100%) - the "sweep" function is helpful for this.
fungi_asv_abun <- data.frame(sweep(fungi_ASVs, 2, colSums(fungi_ASVs), '/')) * 100

#Next aggregate ASV abundances by genus.

fungi_asv_abun$genus <- fungi_taxa_breakdown[rownames(fungi_asv_abun), "genus"]
fungi_asv_abun_genus_sum <- aggregate(. ~ genus, data=fungi_asv_abun, FUN=sum)

#We need to identify rare genera and collapse them into the "Other" category. 

#Reorder based on abundance
fungi_asv_abun_genus_sum <- fungi_asv_abun_genus_sum[order(rowSums(fungi_asv_abun_genus_sum[,2:ncol(fungi_asv_abun_genus_sum)]),decreasing=T),]
#Keep the top 15 but rename everything else starting at the 16th most abundant as Other 
fungi_asv_abun_genus_sum[16:nrow(fungi_asv_abun_genus_sum),1] <- "Other"

#Now melt this table
fungi_asv_abun_relab_genus_sum_melt <- melt(fungi_asv_abun_genus_sum)
#Join in meta_data and change sample names 
fungi_meta <- fungi_meta %>% mutate(variable=as.character(rownames(fungi_meta)))
fungi_asv_abun_relab_genus_sum_melt <-  fungi_asv_abun_relab_genus_sum_melt %>% inner_join(fungi_meta)

#reorder samples based on group
fungi_asv_abun_relab_genus_sum_melt$group <-  gsub("root (cover)", "cover", fungi_asv_abun_relab_genus_sum_melt$group, fixed=T)

fungi_asv_abun_relab_genus_sum_melt$group_num <- as.numeric(as.factor(fungi_asv_abun_relab_genus_sum_melt$group))

fungi_asv_abun_relab_genus_sum_melt$variable <- fct_reorder(fungi_asv_abun_relab_genus_sum_melt$variable, fungi_asv_abun_relab_genus_sum_melt$group_num)

fungi_asv_abun_relab_genus_sum_melt$variable  <- fct_relevel(fungi_asv_abun_relab_genus_sum_melt$variable ,"e98", "e99")

#get custom colour palette
colour_count = length(unique(fungi_asv_abun_relab_genus_sum_melt$genus))
cols <- inlmisc::GetColors(n = 16, scheme = "discrete rainbow")
cols[colour_count] <- "darkgrey"
my_palette <- as.vector(cols[1:16])


# Clean up genus labels:
fungi_asv_abun_relab_genus_sum_melt$genus <- as.character(fungi_asv_abun_relab_genus_sum_melt$genus)
fungi_asv_abun_relab_genus_sum_melt$genus <- gsub("k__Fungi;p__", "", fungi_asv_abun_relab_genus_sum_melt$genus)
fungi_asv_abun_relab_genus_sum_melt$genus <- gsub("k__Unclassified.*$", "Unclassified", fungi_asv_abun_relab_genus_sum_melt$genus)
fungi_asv_abun_relab_genus_sum_melt$genus <- gsub(";c__", "_", fungi_asv_abun_relab_genus_sum_melt$genus)
fungi_asv_abun_relab_genus_sum_melt$genus <- gsub(";o__", "_", fungi_asv_abun_relab_genus_sum_melt$genus)
fungi_asv_abun_relab_genus_sum_melt$genus <- gsub(";f__", "_", fungi_asv_abun_relab_genus_sum_melt$genus)
fungi_asv_abun_relab_genus_sum_melt$genus <- gsub(";g__", "_", fungi_asv_abun_relab_genus_sum_melt$genus)

#reorder genus based on abudance (low to high, across all samples)
fungi_asv_abun_relab_genus_sum_melt$genus <- fct_reorder(fungi_asv_abun_relab_genus_sum_melt$genus, fungi_asv_abun_relab_genus_sum_melt$value, sum)

#Put Other last
fungi_asv_abun_relab_genus_sum_melt$genus <- fct_relevel(fungi_asv_abun_relab_genus_sum_melt$genus, "Other", after = Inf)

fungi_x_col <- rep(NA, length(levels(fungi_asv_abun_relab_genus_sum_melt$variable)))
names(fungi_x_col) <- levels(fungi_asv_abun_relab_genus_sum_melt$variable)
fungi_x_col[fungi_meta[which(fungi_meta$group == "root"), "variable"]] <- "#009E73"
fungi_x_col[fungi_meta[which(fungi_meta$group == "root (cover)"), "variable"]] <- "#E69F00"
fungi_x_col[fungi_meta[which(fungi_meta$group == "soil"), "variable"]] <- "#56B4E9"

#plot
fungi_stacked <- ggplot(fungi_asv_abun_relab_genus_sum_melt, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("Fungi")+
  theme(legend.position="bottom",
        legend.key.size = unit(1, "cm"),
        legend.title=element_blank(),
        axis.text.x = element_text(size=20,
                                   angle = 45,
                                   hjust = 1,
                                   colour=fungi_x_col),
        axis.text.y = element_text(size=20),
        legend.text=element_text(size=17),
        axis.title=element_text(size=20),
        title=element_text(size=20)) +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  scale_fill_manual(values = my_palette)

###Figure 1####
pdf("figure1.pdf", width=8, height=13,family="Arial")
first_row = plot_grid(grobTree(bacteria_sampletype_venn),grobTree(fungi_sampletype_venn), labels = c('A', 'B'), label_size=25)
second_row = plot_grid(bacteria_stacked, labels = c('C'), nrow = 1, label_size=25)
third_row = plot_grid(fungi_stacked, labels = c('D'), nrow = 1, label_size=25)
plot_grid(first_row, second_row, third_row, labels=c('', '', ''), ncol=1)
dev.off()
