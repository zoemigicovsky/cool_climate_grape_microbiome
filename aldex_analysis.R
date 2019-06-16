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
root_depth_aldex_kw <- cbind(rownames(root_depth_aldex_kw), root_depth_aldex_kw)
colnames(root_depth_aldex_kw)[1] <- "bacteria genus"
write.table(root_depth_aldex_kw, "bacteria_root_depth_aldex_kw.txt", sep="\t", quote=F, row.names = F)

rootstock_aldex <- aldex.clr(bacteria_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw <- cbind(rownames(rootstock_aldex_kw), rootstock_aldex_kw)
colnames(rootstock_aldex_kw)[1] <- "bacteria genus"
write.table(rootstock_aldex_kw, "bacteria_rootstock_aldex_kw.txt", sep="\t", quote=F, row.names = F)

#combine for supplemental table
root_depth_aldex_kw <- read_tsv("bacteria_root_depth_aldex_kw.txt")
root_depth_aldex_kw <- root_depth_aldex_kw %>% mutate(factor="root depth") %>% select("bacteria genus", factor, kw.ep:glm.eBH)

rootstock_aldex_kw <- read_tsv("bacteria_rootstock_aldex_kw.txt")
rootstock_aldex_kw <- rootstock_aldex_kw %>% mutate(factor="rootstock") %>% select("bacteria genus", factor, kw.ep:glm.eBH)

aldex_bacteria_results <- bind_rows(root_depth_aldex_kw,rootstock_aldex_kw)
write.table(aldex_bacteria_results, "aldex_bacteria_results.csv", sep=",", quote=F, row.names = F)

#Filter down to only significant
aldex_bacteria_results_sig <- aldex_bacteria_results[aldex_bacteria_results[,"kw.eBH"] <=0.05,]

#In this case there would be none that are significant across root depths
#12 that are significant across rootstocks
write.table(aldex_bacteria_results_sig, "aldex_bacteria_results_sig.csv", sep=",", quote=F, row.names = F)

#Look at distributions for significant genera
#Convert to relative abundances first
bacteria_asv_abun <- data.frame(sweep(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)], 2, colSums(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)]), '/')) * 100
rownames(bacteria_asv_abun) <- bacteria_asv_abun_genus_sum[,"genus"]

aldex_bacteria_results_sig <- as.data.frame(aldex_bacteria_results_sig)
bacteria_asv_abun_rootstock <- bacteria_asv_abun[rownames(bacteria_asv_abun) %in% aldex_bacteria_results_sig[,"bacteria genus"],]
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

#Include unique vine info

#I'm guessing it's the "Other info" R14P10V2 info?
bacteria_asv_abun_rootstock_meta <- bacteria_asv_abun_rootstock_meta %>% mutate(other_info=str_replace(other_info, "NYM Own Root ", ""),other_info=str_replace(other_info, "NYM 3309", ""),other_info=str_replace(other_info, "NYM Rg", ""),other_info=str_replace(other_info, "0-15cm", ""),other_info=str_replace(other_info, "15-30cm", ""),other_info=str_replace(other_info, "30-50cm", ""),other_info=str_replace(other_info, "NYM", ""),other_info=str_trim(other_info))

table(bacteria_asv_abun_rootstock_meta$other_info)
#R14P10V2 R14P10V4 R15P18V1 R15P18V2 R16P19V1 R16P19V3  R17P6V1  R17P6V4 R20P20V2 R20P20V3  R21P5V1  R21P5V3 
#36       36       36       36       24       36       36       36       24       36       36       36  

pdf("rootstock_bacteria_distributions.pdf", width=6.5, height=8,family="Arial")
bacteria_asv_abun_rootstock_meta %>% ggplot(aes(x=factor(rootstock, level=c("Ungrafted", "3309 C", "Riparia Gloire")), y=relative_abun)) + geom_quasirandom(alpha=0.7, stroke=0, size=2) + geom_boxplot(aes(fill=rootstock), alpha=0.3, outlier.alpha = 0) + facet_wrap(~genus_name, scales = "free", ncol=4) +theme_few()+theme(legend.position = "none")+ scale_fill_manual(values = color_palette)+labs(y="Relative Abundance (%)", x="Rootstock")+theme(axis.text = element_text(size=8, colour="black", face="plain"), text=element_text(size=9, face="bold"), axis.text.x = element_text(angle=45, hjust=1, color="black"))
dev.off()

my_palette <- c("#5289C7", "#E8601C", "#4EB265", "#F4A736" ,"#F7F056" ,"#1965B0", "#882E72" ,"#AE76A3", "#CAE0AB" ,"#D1BBD7" , "#DC050C", "#7BAFDE")

pdf("rootstock_bacteria_distributions_vine.pdf", width=7, height=10,family="Arial")
ggplot(bacteria_asv_abun_rootstock_meta, aes(x=factor(rootstock, level=c("Ungrafted", "3309 C", "Riparia Gloire")), y=relative_abun)) + geom_quasirandom(aes(colour=other_info), stroke=0, alpha=0.7, size=3)+ scale_colour_manual(values = my_palette) + geom_boxplot(alpha=0.3, outlier.alpha = 0) + facet_wrap(~genus_name, scales = "free") +theme_few()+theme(legend.position = "bottom")+labs(y="Relative Abundance (%)", x="Rootstock")+theme(axis.text = element_text(size=8, colour="black", face="plain"), text=element_text(size=9, face="bold"), axis.text.x = element_text(angle=45, hjust=1, color="black"))
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
root_depth_aldex_kw <- cbind(rownames(root_depth_aldex_kw), root_depth_aldex_kw)
colnames(root_depth_aldex_kw)[1] <- "fungi genus"
write.table(root_depth_aldex_kw, "fungi_root_depth_aldex_kw.txt", sep="\t", quote=F, row.names = F)

rootstock_aldex <- aldex.clr(fungi_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw <- cbind(rownames(rootstock_aldex_kw), rootstock_aldex_kw)
colnames(rootstock_aldex_kw)[1] <- "fungi genus"
write.table(rootstock_aldex_kw, "fungi_rootstock_aldex_kw.txt", sep="\t", quote=F, row.names = F)

#combine for supplemental table
root_depth_aldex_kw <- read_tsv("fungi_root_depth_aldex_kw.txt")
root_depth_aldex_kw <- root_depth_aldex_kw %>% mutate(factor="root depth") %>% select("fungi genus", factor, kw.ep:glm.eBH)

rootstock_aldex_kw <- read_tsv("fungi_rootstock_aldex_kw.txt")
rootstock_aldex_kw <- rootstock_aldex_kw %>% mutate(factor="rootstock") %>% select("fungi genus", factor, kw.ep:glm.eBH)

aldex_fungi_results <- bind_rows(root_depth_aldex_kw,rootstock_aldex_kw)
write.table(aldex_fungi_results, "aldex_fungi_results.csv", sep=",", quote=F, row.names = F)

#Filter down to only significant
aldex_fungi_results_sig <- aldex_fungi_results[aldex_fungi_results[,"kw.eBH"] <=0.05,]

#In this case there would be none that are significant 
