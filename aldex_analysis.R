#Determine rankings for features for both bacteria and fungi, root depth and rootstock (4 panel figure)

library(ALDEx2)

#The default aldex2 command is for 2 sample groupings but you can use “aldex.kw” which is for 2 or more groups. 
#aldex.kw calculates the expected values of the Kruskal-Wallis test and a glm ANOVA on the data returned by aldex.clr.

source("root_depth_project_functions.R")

###BACTERIA####
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

#Next aggregate ASV abundances by genus.
bacteria_ASVs$genus <- bacteria_taxa_breakdown[rownames(bacteria_ASVs), "genus"]
bacteria_asv_abun_genus_sum <- aggregate(. ~ genus, data=bacteria_ASVs, FUN=sum)

#Convert genus info to rownames 
bacteria_asv_aldex <- bacteria_asv_abun_genus_sum
rownames(bacteria_asv_aldex) <- bacteria_asv_aldex[,1]
bacteria_asv_aldex <- bacteria_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(bacteria_meta)==colnames(bacteria_asv_aldex))

rootstock <- bacteria_meta[,"rootstock"]
root_depth <- bacteria_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(bacteria_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
write.table(root_depth_aldex_kw, "bacteria_root_depth_aldex_kw.txt", sep="\t", quote=F, row.names = T)

rootstock_aldex <- aldex.clr(bacteria_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
write.table(rootstock_aldex_kw, "bacteria_rootstock_aldex_kw.txt", sep="\t", quote=F, row.names = T)

root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#In this case there would be none that are significant across root depths
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#9 that are significant across rootstocks
write.table(rootstock_aldex_kw_sig, "bacteria_rootstock_aldex_kw_sig_genus.txt", sep="\t", quote=F, row.names = T)

#Look at distributions for significant genera
#Convert to relative abundances first
bacteria_asv_abun <- data.frame(sweep(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)], 2, colSums(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)]), '/')) * 100
rownames(bacteria_asv_abun) <- bacteria_asv_abun_genus_sum[,"genus"]

bacteria_asv_abun_rootstock <- bacteria_asv_abun[rownames(bacteria_asv_abun) %in% rownames(rootstock_aldex_hw_sig),]
bacteria_asv_abun_rootstock <- cbind(rownames(bacteria_asv_abun_rootstock), bacteria_asv_abun_rootstock)
colnames(bacteria_asv_abun_rootstock)[1] <- "genus"

library(tidyverse)
library(ggbeeswarm)
library(ggthemes)
library(scales)
library(extrafont)
color_palette <- colorblind_pal()(8)[2:4]
bacteria_asv_abun_rootstock <- bacteria_asv_abun_rootstock %>% gather(key="sample", value="relative_abun", -genus)
#Join with metadata 
bacteria_meta <- cbind(rownames(bacteria_meta), bacteria_meta)
colnames(bacteria_meta)[1] <- "sample"
bacteria_asv_abun_rootstock_meta <- bacteria_asv_abun_rootstock %>% left_join(bacteria_meta)
#get just genus name
bacteria_asv_abun_rootstock_meta <- bacteria_asv_abun_rootstock_meta %>% mutate(genus_name = gsub("D_0__Bacteria;D_1__","",genus))

#Plot distributions for genera that differ significantly 

#rename rootstocks for plotting
bacteria_asv_abun_rootstock_meta <-  bacteria_asv_abun_rootstock_meta %>% mutate(rootstock=str_replace(rootstock, "new_york_muscat", "Ungrafted"), rootstock=str_replace(rootstock, "c3309", "3309 C"),rootstock=str_replace(rootstock, "riparia_gloire", "Riparia Gloire"))
       
pdf("rootstock_bacteria_distributions.pdf", width=6.5, height=8,family="Arial")
bacteria_asv_abun_rootstock_meta %>% ggplot(aes(x=factor(rootstock, level=c("Ungrafted", "3309 C", "Riparia Gloire")), y=relative_abun)) + geom_quasirandom(alpha=0.7, stroke=0, size=2) + geom_boxplot(aes(fill=rootstock), alpha=0.3, outlier.alpha = 0) + facet_wrap(~genus_name, scales = "free") +theme_few()+theme(legend.position = "none")+ scale_fill_manual(values = color_palette)+labs(y="Relative Abundance (%)", x="Rootstock")+theme(axis.text = element_text(size=8, colour="black", face="plain"), text=element_text(size=9, face="bold"), axis.text.x = element_text(angle=45, hjust=1, color="black"))
dev.off()

###FUNGI#### 
fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#looad in grape root only ASVs
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

#Next aggregate ASV abundances by genus.
fungi_ASVs$genus <- fungi_taxa_breakdown[rownames(fungi_ASVs), "genus"]
fungi_asv_abun_genus_sum <- aggregate(. ~ genus, data=fungi_ASVs, FUN=sum)

#Convert genus info to rownames 
fungi_asv_aldex <- fungi_asv_abun_genus_sum
rownames(fungi_asv_aldex) <- fungi_asv_aldex[,1]
fungi_asv_aldex <- fungi_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(fungi_meta)==colnames(fungi_asv_aldex))

rootstock <- fungi_meta[,"rootstock"]
root_depth <- fungi_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(fungi_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
write.table(root_depth_aldex_kw, "fungi_root_depth_aldex_kw.txt", sep="\t", quote=F, row.names = T)

rootstock_aldex <- aldex.clr(fungi_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
write.table(rootstock_aldex_kw, "fungi_rootstock_aldex_kw.txt", sep="\t", quote=F, row.names = T)

root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#In this case there would be none that are significant across root depths
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#Based on glm.eBH there are no significant genera for fungi 
