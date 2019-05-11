library(readxl)
library(tidyverse)

soil_dat <- read_excel("(2018) NYM rootstock depth trial all data updated.xlsx", skip = 4) %>% select("...3", pH, "Organic Matter":"Mg...14", "Na...19":H)

soil_dat <- soil_dat[-1,]

colnames(soil_dat) <- c("id", "ph", "organic_matter", "p2o5", "k2o",  "ca", "mg", "na", "s", "al", "b", "cu", "fe", "mn", "zn", "cec", "k_sat", "ca_sat", "mg_sat", "na_sat", "h_sat")

soil_dat <- soil_dat %>% filter( id != "NA")

soil_dat <- soil_dat %>% mutate(b=str_replace(b, "< 0.50","0.25"))

write.table(soil_dat, "soil_dat.tsv", sep="\t", row.names=F, quote=F, col.names=T)

soil_dat <- read_tsv("soil_dat.tsv")

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
spearman_results_diversity <- spearman_results %>% mutate(p.value=p.value*20) %>% mutate(p.value=ifelse(p.value>1, 1, p.value))

write.table(spearman_results_diversity, "bacteria_soil_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Now do the same thing but with the main genera 

#FUNGI ANALYSIS 

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
spearman_results_diversity <- spearman_results %>% mutate(p.value=p.value*20) %>% mutate(p.value=ifelse(p.value>1, 1, p.value))

write.table(spearman_results_diversity, "fungi_soil_diversity_spearman_results.csv", sep=",", col.names = T, row.names=F, quote=F)

#Now do the same thing but with the main genera 

