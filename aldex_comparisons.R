library(ALDEx2)

#ASVS

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

#Need a character vector of group labels
table(rownames(bacteria_meta)==colnames(bacteria_ASVs))

rootstock <- bacteria_meta[,"rootstock"]
root_depth <- bacteria_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(bacteria_ASVs, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
write.table(root_depth_aldex_kw, "bacteria_root_depth_aldex_kw_asv.txt", sep="\t", quote=F, row.names = T)

root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#1 sig

rootstock_aldex <- aldex.clr(bacteria_ASVs, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
write.table(rootstock_aldex_kw, "bacteria_rootstock_aldex_kw_asv.txt", sep="\t", quote=F, row.names = T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#2 that are significant across rootstocks

###FUNGI#### 

fungi_meta <- read.table("fungi/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#load in grape root only ASVs
fungi_ASVs <- read.table("fungi/dada2_output_exported_grape/feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

#Keep only meta data for grape roots
fungi_meta <- fungi_meta[colnames(fungi_ASVs),]

table(colnames(fungi_ASVs)== rownames(fungi_meta))
#they match

#Need a character vector of group labels
table(rownames(fungi_meta)==colnames(fungi_ASVs))

rootstock <- fungi_meta[,"rootstock"]
root_depth <- fungi_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(fungi_ASVs, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
write.table(root_depth_aldex_kw, "fungi_root_depth_aldex_kw_asv.txt", sep="\t", quote=F, row.names = T)

root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#0 sig

rootstock_aldex <- aldex.clr(fungi_ASVs, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
write.table(rootstock_aldex_kw, "fungi_rootstock_aldex_kw_asv.txt", sep="\t", quote=F, row.names = T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#0 sig

#Redo analysis with family-level taxonomy 

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

#Next aggregate ASV abundances by family
bacteria_ASVs$family <- bacteria_taxa_breakdown[rownames(bacteria_ASVs), "family"]
bacteria_asv_abun_family_sum <- aggregate(. ~ family, data=bacteria_ASVs, FUN=sum)

#Convert family info to rownames 
bacteria_asv_aldex <- bacteria_asv_abun_family_sum
rownames(bacteria_asv_aldex) <- bacteria_asv_aldex[,1]
bacteria_asv_aldex <- bacteria_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(bacteria_meta)==colnames(bacteria_asv_aldex))

rootstock <- bacteria_meta[,"rootstock"]
root_depth <- bacteria_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(bacteria_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#Nothing significant 

rootstock_aldex <- aldex.clr(bacteria_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#7 that are significant across rootstocks
write.table(rootstock_aldex_kw_sig, "bacteria_rootstock_aldex_kw_sig_family.txt", sep="\t", quote=F, row.names = T)

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
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

#Next aggregate ASV abundances by family.
fungi_ASVs$family <- fungi_taxa_breakdown[rownames(fungi_ASVs), "family"]
fungi_asv_abun_family_sum <- aggregate(. ~ family, data=fungi_ASVs, FUN=sum)

#Convert family info to rownames 
fungi_asv_aldex <- fungi_asv_abun_family_sum
rownames(fungi_asv_aldex) <- fungi_asv_aldex[,1]
fungi_asv_aldex <- fungi_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(fungi_meta)==colnames(fungi_asv_aldex))

rootstock <- fungi_meta[,"rootstock"]
root_depth <- fungi_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(fungi_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#1 sig

rootstock_aldex <- aldex.clr(fungi_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#None

###ORDER

library(ALDEx2)

#Redo analysis with order-level taxonomy 

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

#Next aggregate ASV abundances by order
bacteria_ASVs$order <- bacteria_taxa_breakdown[rownames(bacteria_ASVs), "order"]
bacteria_asv_abun_order_sum <- aggregate(. ~ order, data=bacteria_ASVs, FUN=sum)

#Convert order info to rownames 
bacteria_asv_aldex <- bacteria_asv_abun_order_sum
rownames(bacteria_asv_aldex) <- bacteria_asv_aldex[,1]
bacteria_asv_aldex <- bacteria_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(bacteria_meta)==colnames(bacteria_asv_aldex))

rootstock <- bacteria_meta[,"rootstock"]
root_depth <- bacteria_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(bacteria_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#Nothing significant 

rootstock_aldex <- aldex.clr(bacteria_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#4 that are significant across rootstocks
write.table(rootstock_aldex_kw_sig, "bacteria_rootstock_aldex_kw_sig_order.txt", sep="\t", quote=F, row.names = T)

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
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

#Next aggregate ASV abundances by order.
fungi_ASVs$order <- fungi_taxa_breakdown[rownames(fungi_ASVs), "order"]
fungi_asv_abun_order_sum <- aggregate(. ~ order, data=fungi_ASVs, FUN=sum)

#Convert order info to rownames 
fungi_asv_aldex <- fungi_asv_abun_order_sum
rownames(fungi_asv_aldex) <- fungi_asv_aldex[,1]
fungi_asv_aldex <- fungi_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(fungi_meta)==colnames(fungi_asv_aldex))

rootstock <- fungi_meta[,"rootstock"]
root_depth <- fungi_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(fungi_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#1 sig

rootstock_aldex <- aldex.clr(fungi_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#None

###CLASS

library(ALDEx2)

#Redo analysis with class-level taxonomy 

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

#Next aggregate ASV abundances by class
bacteria_ASVs$class <- bacteria_taxa_breakdown[rownames(bacteria_ASVs), "class"]
bacteria_asv_abun_class_sum <- aggregate(. ~ class, data=bacteria_ASVs, FUN=sum)

#Convert class info to rownames 
bacteria_asv_aldex <- bacteria_asv_abun_class_sum
rownames(bacteria_asv_aldex) <- bacteria_asv_aldex[,1]
bacteria_asv_aldex <- bacteria_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(bacteria_meta)==colnames(bacteria_asv_aldex))

rootstock <- bacteria_meta[,"rootstock"]
root_depth <- bacteria_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(bacteria_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#Nothing significant 

rootstock_aldex <- aldex.clr(bacteria_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#0 that are significant across rootstocks

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
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

#Next aggregate ASV abundances by class.
fungi_ASVs$class <- fungi_taxa_breakdown[rownames(fungi_ASVs), "class"]
fungi_asv_abun_class_sum <- aggregate(. ~ class, data=fungi_ASVs, FUN=sum)

#Convert class info to rownames 
fungi_asv_aldex <- fungi_asv_abun_class_sum
rownames(fungi_asv_aldex) <- fungi_asv_aldex[,1]
fungi_asv_aldex <- fungi_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(fungi_meta)==colnames(fungi_asv_aldex))

rootstock <- fungi_meta[,"rootstock"]
root_depth <- fungi_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(fungi_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#1 significant 

rootstock_aldex <- aldex.clr(fungi_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#None

#Redo analysis with phylum-level taxonomy 

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

#Next aggregate ASV abundances by phylum
bacteria_ASVs$phylum <- bacteria_taxa_breakdown[rownames(bacteria_ASVs), "phylum"]
bacteria_asv_abun_phylum_sum <- aggregate(. ~ phylum, data=bacteria_ASVs, FUN=sum)

#Convert phylum info to rownames 
bacteria_asv_aldex <- bacteria_asv_abun_phylum_sum
rownames(bacteria_asv_aldex) <- bacteria_asv_aldex[,1]
bacteria_asv_aldex <- bacteria_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(bacteria_meta)==colnames(bacteria_asv_aldex))

rootstock <- bacteria_meta[,"rootstock"]
root_depth <- bacteria_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(bacteria_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#Nothing significant 

rootstock_aldex <- aldex.clr(bacteria_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#1 that are significant across rootstocks

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
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

#Next aggregate ASV abundances by phylum.
fungi_ASVs$phylum <- fungi_taxa_breakdown[rownames(fungi_ASVs), "phylum"]
fungi_asv_abun_phylum_sum <- aggregate(. ~ phylum, data=fungi_ASVs, FUN=sum)

#Convert phylum info to rownames 
fungi_asv_aldex <- fungi_asv_abun_phylum_sum
rownames(fungi_asv_aldex) <- fungi_asv_aldex[,1]
fungi_asv_aldex <- fungi_asv_aldex[,-1]

#Need a character vector of group labels
table(rownames(fungi_meta)==colnames(fungi_asv_aldex))

rootstock <- fungi_meta[,"rootstock"]
root_depth <- fungi_meta[,"root_depth"]

root_depth_aldex <- aldex.clr(fungi_asv_aldex, root_depth,mc.samples = 128, verbose=T)
root_depth_aldex_kw <- aldex.kw(root_depth_aldex, verbose=T)
root_depth_aldex_kw_sig <- root_depth_aldex_kw[root_depth_aldex_kw[,"glm.eBH"] <=0.05,]
#None

rootstock_aldex <- aldex.clr(fungi_asv_aldex, rootstock,mc.samples = 128, verbose=T)
rootstock_aldex_kw <- aldex.kw(rootstock_aldex, verbose=T)
rootstock_aldex_kw_sig <- rootstock_aldex_kw[rootstock_aldex_kw[,"glm.eBH"] <=0.05,]
#None

