library(tidyverse)
library(broom)
#Bacteria correlation with root phenotypes 

root_phenos <- read_tsv("root_phenotypes.txt") %>% filter(!is.na(sample_id))

#Alpha diversity - richness 

alpha_diversity <- read_tsv("bacteria/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_bacteria <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_bacteria <- root_phenos_bacteria %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(observed_otus, trunk_diameter:total_root_mass)

#Run correlations with soil_fungi$observed_otus
root_cor <- apply(root_phenos_bacteria[, -1], 2, cor.test, root_phenos_bacteria$observed_otus, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(root_cor)){
  spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests 
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "bacteria_root_pheno_richness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Evenness 
alpha_diversity <- read_tsv("bacteria/diversity_grape/evenness_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_bacteria <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_bacteria <- root_phenos_bacteria %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(pielou_e, trunk_diameter:total_root_mass)

#Run correlations with soil_fungi$observed_otus
root_cor <- apply(root_phenos_bacteria[, -1], 2, cor.test, root_phenos_bacteria$pielou_e, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(root_cor)){
  spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests 
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "bacteria_root_pheno_evenness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Faiths 
alpha_diversity <- read_tsv("bacteria/diversity_grape/faith_pd_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_bacteria <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_bacteria <- root_phenos_bacteria %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(faith_pd, trunk_diameter:total_root_mass)

#Run correlations with soil_fungi$observed_otus
root_cor <- apply(root_phenos_bacteria[, -1], 2, cor.test, root_phenos_bacteria$faith_pd, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(root_cor)){
  spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests 
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "bacteria_root_pheno_faith_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)


###now do the same thing with fungi

root_phenos <- read_tsv("root_phenotypes.txt") %>% filter(!is.na(sample_id))

#Alpha diversity - richness 

alpha_diversity <- read_tsv("fungi/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_fungi <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_fungi <- root_phenos_fungi %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(observed_otus, trunk_diameter:total_root_mass)

#Run correlations with soil_fungi$observed_otus
root_cor <- apply(root_phenos_fungi[, -1], 2, cor.test, root_phenos_fungi$observed_otus, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(root_cor)){
  spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests 
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "fungi_root_pheno_richness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Evenness 
alpha_diversity <- read_tsv("fungi/diversity_grape/evenness_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_fungi <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_fungi <- root_phenos_fungi %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(pielou_e, trunk_diameter:total_root_mass)

#Run correlations with soil_fungi$observed_otus
root_cor <- apply(root_phenos_fungi[, -1], 2, cor.test, root_phenos_fungi$pielou_e, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(root_cor)){
  spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests 
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "fungi_root_pheno_evenness_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Faiths 
alpha_diversity <- read_tsv("fungi/diversity_grape/faith_pd_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_fungi <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_fungi <- root_phenos_fungi %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(faith_pd, trunk_diameter:total_root_mass)

#Run correlations with soil_fungi$observed_otus
root_cor <- apply(root_phenos_fungi[, -1], 2, cor.test, root_phenos_fungi$faith_pd, method="spearman")

spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )

for(i in 1:length(root_cor)){
  spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
  spearman_results <- rbind(spearman_results, spear_test)
}
spearman_results <- spearman_results[-1,]
spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)

#Multiple p-value by number of tests 
spearman_results$p.value <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))

write.table(spearman_results, "fungi_root_pheno_faith_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

