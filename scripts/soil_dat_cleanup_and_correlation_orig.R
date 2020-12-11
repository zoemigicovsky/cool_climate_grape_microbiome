library(readxl)
library(tidyverse)
library(broom)

setwd("data/")

soil_dat <- read_excel("../(2018) NYM rootstock depth trial all data updated.xlsx", skip = 4) %>% select("...3", pH, "Organic Matter":"Mg...14", "Na...19":H)

soil_dat <- soil_dat[-1,]

colnames(soil_dat) <- c("id", "ph", "organic_matter", "p2o5", "k2o",  "ca", "mg", "na", "s", "al", "b", "cu", "fe", "mn", "zn", "cec", "k_sat", "ca_sat", "mg_sat", "na_sat", "h_sat")

soil_dat <- soil_dat %>% filter( id != "NA")

soil_dat <- soil_dat %>% mutate(b=str_replace(b, "< 0.50","0.25"))

write.table(soil_dat, "soil_dat.tsv", sep="\t", row.names=F, quote=F, col.names=T)

#BACTERIA ANALYSIS 

#Alpha diversity

alpha_diversity <- read_tsv("bacteria/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

soil_bacteria <- soil_dat %>% rename(sample_id=id) %>% inner_join(alpha_diversity) %>% select(observed_otus, ph:h_sat)

#Run correlations with soil_bacteria$observed_otus
soil_cor <- apply(soil_bacteria[, -1], 2, cor.test, soil_bacteria$observed_otus, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(soil_cor)){
  spear_test <- tidy(soil_cor[[i]]) %>% mutate(pheno=names(soil_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "bacteria_soil_richness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Faith
alpha_diversity <- read_tsv("bacteria/diversity_grape/faith_pd_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

soil_bacteria <- soil_dat %>% rename(sample_id=id) %>% inner_join(alpha_diversity) %>% select(faith_pd, ph:h_sat)

#Run correlations with soil_bacteria$faith_pd
soil_cor <- apply(soil_bacteria[, -1], 2, cor.test, soil_bacteria$faith_pd, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(soil_cor)){
  spear_test <- tidy(soil_cor[[i]]) %>% mutate(pheno=names(soil_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "bacteria_soil_faith_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Evenness

alpha_diversity <- read_tsv("bacteria/diversity_grape/evenness_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

soil_bacteria <- soil_dat %>% rename(sample_id=id) %>% inner_join(alpha_diversity) %>% select(pielou_e, ph:h_sat)

#Run correlations with soil_bacteria$faith_pd
soil_cor <- apply(soil_bacteria[, -1], 2, cor.test, soil_bacteria$pielou_e, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(soil_cor)){
  spear_test <- tidy(soil_cor[[i]]) %>% mutate(pheno=names(soil_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "bacteria_soil_evenness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Combine for a supplemental table 

bacteria_evenness <- read_csv("bacteria_soil_evenness_diversity_spearman_results.csv") 
bacteria_faith <- read_csv("bacteria_soil_faith_diversity_spearman_results.csv") 
bacteria_richness <- read_csv("bacteria_soil_richness_diversity_spearman_results.csv") 

bacteria_evenness <- bacteria_evenness %>% mutate(diversity_metric="evenness")
bacteria_faith <- bacteria_faith %>% mutate(diversity_metric="Faith's Phylogenetic Diversity")
bacteria_richness <- bacteria_richness %>% mutate(diversity_metric="richness")

bacteria_diversity <- bind_rows(bacteria_evenness,bacteria_faith,bacteria_richness) %>% select(diversity_metric, pheno, estimate, p.value)

write.table(bacteria_diversity, "bacteria_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)


#Now do the same thing but with the main genera 

#load in raw feature data 
library(reshape2)
library(tidyverse)
library(RColorBrewer)
source("root_depth_project_functions.R")

bacteria_meta <- read.table("bacteria/root_depth_bacteria_metadata_grape.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

bacteria_ASVs <- read.table("bacteria/dada2_output_exported_grape//feature-table_w_tax.txt",
                            header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

bacteria_ASVs <- bacteria_ASVs[, -which(colnames(bacteria_ASVs) == "taxonomy")]

bacteria_taxa <- read.table("ASV_tables/bacteria/taxonomy.tsv",
                            header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
bacteria_taxa_breakdown <- qiime2_taxa_breakdown(bacteria_taxa)

#Convert to relative abundance (i.e. sum to 100%) - the "sweep" function is helpful for this.
bacteria_asv_abun <- data.frame(sweep(bacteria_ASVs, 2, colSums(bacteria_ASVs), '/')) * 100

#Next aggregate ASV abundances by genus.

bacteria_asv_abun$genus <- bacteria_taxa_breakdown[rownames(bacteria_asv_abun), "genus"]
bacteria_asv_abun_genus_sum <- aggregate(. ~ genus, data=bacteria_asv_abun, FUN=sum)

#Get the list and data for the 10 most abundant grape genera are

#Reorder based on abundance
bacteria_asv_abun_genus_sum <- bacteria_asv_abun_genus_sum[order(rowSums(bacteria_asv_abun_genus_sum[,2:ncol(bacteria_asv_abun_genus_sum)]),decreasing=T),]

#Keep the top 10

bacteria_asv_abun_genus_sum_top <- bacteria_asv_abun_genus_sum[1:10,]
bacteria_asv_abun_genus_sum_top <- as_tibble(cbind(nms = names(bacteria_asv_abun_genus_sum_top), t(bacteria_asv_abun_genus_sum_top)))
colnames(bacteria_asv_abun_genus_sum_top) <- bacteria_asv_abun_genus_sum_top[1,]
bacteria_asv_abun_genus_sum_top <- bacteria_asv_abun_genus_sum_top[-1,]
colnames(bacteria_asv_abun_genus_sum_top)[1] <- "id"
bacteria_asv_abun_genus_sum_top <- bacteria_asv_abun_genus_sum_top %>% arrange(id)

#Correlate with soil

soil_dat <- read_tsv("soil_dat.tsv")

#Reduce phenotype data down to samples with diversity metrics

soil_bacteria <- soil_dat %>% semi_join(bacteria_asv_abun_genus_sum_top) %>% arrange(id) 

table(soil_bacteria[,1] == bacteria_asv_abun_genus_sum_top[,1])

soil_bacteria <- soil_bacteria %>% select(-id)
bacteria_asv_abun_genus_sum_top <- bacteria_asv_abun_genus_sum_top %>% select(-id)

#convert abundances to numberic
bacteria_asv_abun_genus_sum_top <- sapply(bacteria_asv_abun_genus_sum_top, as.numeric)

#Run correlations
library(tidyverse)
library(broom)
spearman_results <- data.frame(matrix(ncol = 7, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno", "bacteria"  )

for(i in 1:ncol(bacteria_asv_abun_genus_sum_top)){
  soil_cor <- apply(soil_bacteria, 2, cor.test, bacteria_asv_abun_genus_sum_top[,i], method="spearman")
  for(j in 1:length(soil_cor)){
    spear_test <- tidy(soil_cor[[j]]) %>% mutate(pheno=names(soil_cor)[j]) %>% mutate(bacteria=colnames( bacteria_asv_abun_genus_sum_top)[i])
    spearman_results <- rbind(spearman_results, spear_test)
  }
}

spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)  

#Multiple p-value by number of tests

spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))
spearman_results <- spearman_results %>% select(bacteria_genus=bacteria, pheno, estimate, p.value)
  
write.table(spearman_results, "bacteria_genus_soil_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#FUNGI ANALYSIS 

library(readxl)
library(tidyverse)
library(broom)
library(broom)

soil_dat <- read_tsv("soil_dat.tsv")
#Alpha diversity

alpha_diversity <- read_tsv("fungi/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

soil_fungi <- soil_dat %>% rename(sample_id=id) %>% inner_join(alpha_diversity) %>% select(observed_otus, ph:h_sat)

#Run correlations with soil_fungi$observed_otus
soil_cor <- apply(soil_fungi[, -1], 2, cor.test, soil_fungi$observed_otus, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(soil_cor)){
  spear_test <- tidy(soil_cor[[i]]) %>% mutate(pheno=names(soil_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "fungi_soil_richness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Faith
alpha_diversity <- read_tsv("fungi/diversity_grape/faith_pd_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

soil_fungi <- soil_dat %>% rename(sample_id=id) %>% inner_join(alpha_diversity) %>% select(faith_pd, ph:h_sat)

#Run correlations with soil_fungi$faith_pd
soil_cor <- apply(soil_fungi[, -1], 2, cor.test, soil_fungi$faith_pd, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(soil_cor)){
  spear_test <- tidy(soil_cor[[i]]) %>% mutate(pheno=names(soil_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "fungi_soil_faith_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Evenness

alpha_diversity <- read_tsv("fungi/diversity_grape/evenness_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

soil_fungi <- soil_dat %>% rename(sample_id=id) %>% inner_join(alpha_diversity) %>% select(pielou_e, ph:h_sat)

#Run correlations with soil_fungi$faith_pd
soil_cor <- apply(soil_fungi[, -1], 2, cor.test, soil_fungi$pielou_e, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(soil_cor)){
  spear_test <- tidy(soil_cor[[i]]) %>% mutate(pheno=names(soil_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "fungi_soil_evenness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Combine for a supplemental table 

fungi_evenness <- read_csv("fungi_soil_evenness_diversity_spearman_results.csv") 
fungi_faith <- read_csv("fungi_soil_faith_diversity_spearman_results.csv") 
fungi_richness <- read_csv("fungi_soil_richness_diversity_spearman_results.csv") 

fungi_evenness <- fungi_evenness %>% mutate(diversity_metric="evenness")
fungi_faith <- fungi_faith %>% mutate(diversity_metric="Faith's Phylogenetic Diversity")
fungi_richness <- fungi_richness %>% mutate(diversity_metric="richness")

fungi_diversity <- bind_rows(fungi_evenness,fungi_faith,fungi_richness) %>% select(diversity_metric, pheno, estimate, p.value)

write.table(fungi_diversity, "fungi_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)


#Now do the same thing but with the main genera 

#load in raw feature data 
source("root_depth_project_functions.R")

fungi_meta <- read.table("metadata/root_depth_fungi_metadata.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

fungi_ASVs <- read.table("fungi/dada2_output_exported_grape//feature-table_w_tax.txt",
                         header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

fungi_ASVs <- fungi_ASVs[, -which(colnames(fungi_ASVs) == "taxonomy")]

fungi_taxa <- read.table("ASV_tables/fungi/taxonomy.tsv",
                         header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors = FALSE)

# Note that the below function to get a dataframe of taxa labels at different levels was sourced from root_depth_project_functions.R.
fungi_taxa_breakdown <- UNITE_qiime2_taxa_breakdown(fungi_taxa)

#Convert to relative abundance (i.e. sum to 100%) - the "sweep" function is helpful for this.
fungi_asv_abun <- data.frame(sweep(fungi_ASVs, 2, colSums(fungi_ASVs), '/')) * 100

#Next aggregate ASV abundances by genus.

fungi_asv_abun$genus <- fungi_taxa_breakdown[rownames(fungi_asv_abun), "genus"]
fungi_asv_abun_genus_sum <- aggregate(. ~ genus, data=fungi_asv_abun, FUN=sum)

#Get the list and data for the 10 most abundant grape genera are

#Reorder based on abundance
fungi_asv_abun_genus_sum <- fungi_asv_abun_genus_sum[order(rowSums(fungi_asv_abun_genus_sum[,2:ncol(fungi_asv_abun_genus_sum)]),decreasing=T),]

#Keep the top 10

fungi_asv_abun_genus_sum_top <- fungi_asv_abun_genus_sum[1:10,]
fungi_asv_abun_genus_sum_top <- as_tibble(cbind(nms = names(fungi_asv_abun_genus_sum_top), t(fungi_asv_abun_genus_sum_top)))
colnames(fungi_asv_abun_genus_sum_top) <- fungi_asv_abun_genus_sum_top[1,]
fungi_asv_abun_genus_sum_top <- fungi_asv_abun_genus_sum_top[-1,]
colnames(fungi_asv_abun_genus_sum_top)[1] <- "id"
fungi_asv_abun_genus_sum_top <- fungi_asv_abun_genus_sum_top %>% arrange(id)

#Correlate with soil

soil_dat <- read_tsv("soil_dat.tsv")

#Reduce phenotype data down to samples with diversity metrics

soil_fungi <- soil_dat %>% semi_join(fungi_asv_abun_genus_sum_top) %>% arrange(id) 

table(soil_fungi[,1] == fungi_asv_abun_genus_sum_top[,1])

soil_fungi <- soil_fungi %>% select(-id)
fungi_asv_abun_genus_sum_top <- fungi_asv_abun_genus_sum_top %>% select(-id)

#convert abundances to numberic
fungi_asv_abun_genus_sum_top <- sapply(fungi_asv_abun_genus_sum_top, as.numeric)

#Run correlations
spearman_results <- data.frame(matrix(ncol = 7, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno", "fungi"  )

for(i in 1:ncol(fungi_asv_abun_genus_sum_top)){
  soil_cor <- apply(soil_fungi, 2, cor.test, fungi_asv_abun_genus_sum_top[,i], method="spearman")
  for(j in 1:length(soil_cor)){
    spear_test <- tidy(soil_cor[[j]]) %>% mutate(pheno=names(soil_cor)[j]) %>% mutate(fungi=colnames( fungi_asv_abun_genus_sum_top)[i])
    spearman_results <- rbind(spearman_results, spear_test)
  }
}

spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)  

#Multiple p-value by number of tests

spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

spearman_results %>% filter(p.value <= 0.05)
#    estimate    p.value  pheno                                                                                fungi
#  0.7366645 0.02119396    na k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pyronemataceae;g__Pseudaleuria
# -0.7103334 0.02993484    fe k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pyronemataceae;g__Pseudaleuria
#  0.7509783 0.02119396    mn k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pyronemataceae;g__Pseudaleuria

spearman_results <- spearman_results %>% select(fungi_genus=fungi, pheno, estimate, p.value)

write.table(spearman_results, "fungi_genus_soil_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

