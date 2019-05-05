#Determine rankings for features for both bacteria and fungi, root depth and rootstock (4 panel figure)

library(ALDEx2)

#The default aldex2 command is for 2 sample groupings but you can use “aldex.kw” which is for 2 or more groups. 

#aldex.kw calculates the expected values of the Kruskal-Wallis test and a glm ANOVA on the data returned by aldex.clr.

#First I need to run aldex.clr then

#Getting started with ALDEx2 is easy. All you need is a matrix (with rows as variables and columns as samples) and a character vector of group labels. Finally, use the denom argument to choose a set of variables to use as the reference for the analysis. 
source("root_depth_project_functions.R")

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#looad in grape root only ASVs
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

#Convert to relative abundance (i.e. sum to 100%) - the "sweep" function is helpful for this.
#bacteria_asv_abun <- data.frame(sweep(bacteria_ASVs, 2, colSums(bacteria_ASVs), '/')) * 100
#I can't get it to run with relative abundances because they have to be integers, try with raw counts

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
rootstock_aldex_hw <- aldex.kw(rootstock_aldex, verbose=T)
write.table(rootstock_aldex_hw, "bacteria_rootstock_aldex_kw.txt", sep="\t", quote=F, row.names = T)

#Returns a data.frame with the following information:
#kw.ep	
#a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#kw.eBH	
#a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#glm.ep	
#a vector containing the expected p-value of the glm ANOVA for each feature
#glm.eBH	
#a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature

root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#In this case there would be none that are significant across root depths
rootstock_aldex_hw_sig <- rootstock_aldex_hw[rootstock_aldex_hw[,"glm.eBH"] <=0.05,]
#9 that are significant across rootstocks
write.table(rootstock_aldex_hw_sig, "bacteria_rootstock_aldex_kw_sig.txt", sep="\t", quote=F, row.names = T)

#Look at distributions for significant genera
#Convert to relative abundances first? 
bacteria_asv_abun <- data.frame(sweep(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)], 2, colSums(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)]), '/')) * 100
rownames(bacteria_asv_abun) <- bacteria_asv_abun_genus_sum[,"genus"]

bacteria_asv_abun_rootstock <- bacteria_asv_abun[rownames(bacteria_asv_abun) %in% rownames(rootstock_aldex_hw_sig),]
bacteria_asv_abun_rootstock <- cbind(rownames(bacteria_asv_abun_rootstock), bacteria_asv_abun_rootstock)
colnames(bacteria_asv_abun_rootstock)[1] <- "genus"

library(tidyverse)
bacteria_asv_abun_rootstock <- bacteria_asv_abun_rootstock %>% gather(key="sample", value="relative_abun", -genus)
#Join with metadata 
bacteria_meta <- cbind(rownames(bacteria_meta), bacteria_meta)
colnames(bacteria_meta)[1] <- "sample"
bacteria_asv_abun_rootstock_meta <- bacteria_asv_abun_rootstock %>% left_join(bacteria_meta)
#get just genus name
bacteria_asv_abun_rootstock_meta <- bacteria_asv_abun_rootstock_meta %>% mutate(genus_name = gsub("D_0__Bacteria;D_1__","",genus))

bacteria_asv_abun_rootstock_meta %>% ggplot(aes(x=rootstock, y=relative_abun, fill=rootstock)) + geom_boxplot() + theme_bw()+theme(legend.position = "none") + facet_wrap(~genus_name, scales = "free")
