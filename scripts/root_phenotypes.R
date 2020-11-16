<<<<<<< HEAD
rm(list=ls(all.names=TRUE))

#Bacteria correlation with root phenotypes 
=======
library(tidyverse)
library(broom)
#Bacteria correlation with root phenotypes 

root_phenos <- read_tsv("data/root_phenotypes.txt") %>% filter(!is.na(sample_id))

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


#Combine for a supplemental table 

bacteria_evenness <- read_csv("bacteria_root_pheno_evenness_diversity_spearman_results.csv") 
bacteria_faith <- read_csv("bacteria_root_pheno_faith_diversity_spearman_results.csv") 
bacteria_richness <- read_csv("bacteria_root_pheno_richness_diversity_spearman_results.csv") 

bacteria_evenness <- bacteria_evenness %>% mutate(diversity_metric="evenness")
bacteria_faith <- bacteria_faith %>% mutate(diversity_metric="Faith's Phylogenetic Diversity")
bacteria_richness <- bacteria_richness %>% mutate(diversity_metric="richness")

bacteria_diversity <- bind_rows(bacteria_evenness,bacteria_faith,bacteria_richness) %>% select(diversity_metric, pheno, estimate, p.value)
>>>>>>> 42510f1db92e8943c212f742b59b31018e899a92

library(tidyverse)
library(broom)

setwd("/home/gavin/github_repos/root_depth/data/")

root_phenos <- read_tsv("data/root_phenotypes.txt") %>% filter(!is.na(sample_id))

# Read in alpha-div restricted to grape root samples only.
bacteria_richness <- read_tsv("diversity_files/bacteria/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)
bacteria_evenness <- read_tsv("diversity_files/bacteria/diversity_grape/evenness_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)
bacteria_phylo_d <- read_tsv("diversity_files/bacteria/diversity_grape/faith_pd_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)
bacteria_shannon <- read_tsv("diversity_files/bacteria/diversity_grape/shannon_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

fungi_richness <- read_tsv("diversity_files/fungi/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)
fungi_evenness <- read_tsv("diversity_files/fungi/diversity_grape/evenness_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)
fungi_phylo_d <- read_tsv("diversity_files/fungi/diversity_grape/faith_pd_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)
fungi_shannon <- read_tsv("diversity_files/fungi/diversity_grape/shannon_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)


root_cor_alpha_and_pheno <- function(root_phenos_in, alpha_div_in, metric_name, kingdom) {

  orig_metric <- metric_name
  
  if(metric_name == "richness") { metric_name <- "observed_otus" }
  if(metric_name == "evenness") { metric_name <- "pielou_e" }
  if(metric_name == "phylo_d") { metric_name <- "faith_pd" }
  
  root_phenos_subset <- root_phenos_in %>% inner_join(alpha_div_in) %>% rename(trunk_diameter='trunk_diameter_at ground') %>% select(metric_name, trunk_diameter:total_root_mass)
  
  root_cor <- apply(root_phenos_subset[, -1], 2, cor.test, root_phenos_subset %>% pull(metric_name), method="spearman")
  
  spearman_results <- data.frame(matrix(ncol = 6, nrow = 1))
  colnames(spearman_results) <- c("estimate",  "statistic" ,  "p.value"  ,   "method"   ,   "alternative", "pheno"  )
  for(i in 1:length(root_cor)){
    spear_test <- tidy(root_cor[[i]]) %>% mutate(pheno=names(root_cor)[i])
    spearman_results <- rbind(spearman_results, spear_test)
  }
  
  spearman_results <- spearman_results[-1,]
  spearman_results <- spearman_results %>% select(-method, -alternative, -statistic)
  
  spearman_results$diversity_metric <- orig_metric
  spearman_results$kingdom <- kingdom
  spearman_results$BH <- p.adjust(spearman_results$p.value, method = "BH", n = length(spearman_results$p.value))
  
  colnames(spearman_results)[which(colnames(spearman_results) == "estimate")] <- "rho"
  
  return(spearman_results[, c("kingdom", "diversity_metric", "pheno", "rho", "p.value", "BH")])
  
}
<<<<<<< HEAD
=======
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

#Combine for a supplemental table 

fungi_evenness <- read_csv("fungi_root_pheno_evenness_diversity_spearman_results.csv") 
fungi_faith <- read_csv("fungi_root_pheno_faith_diversity_spearman_results.csv") 
fungi_richness <- read_csv("fungi_root_pheno_richness_diversity_spearman_results.csv") 

fungi_evenness <- fungi_evenness %>% mutate(diversity_metric="evenness")
fungi_faith <- fungi_faith %>% mutate(diversity_metric="Faith's Phylogenetic Diversity")
fungi_richness <- fungi_richness %>% mutate(diversity_metric="richness")
>>>>>>> 42510f1db92e8943c212f742b59b31018e899a92

input_params <- list()
input_params[["bacteria_richness"]] <- list(alpha_div_in=bacteria_richness, metric_name="richness", kingdom="bacteria")
input_params[["bacteria_evenness"]] <- list(alpha_div_in=bacteria_evenness, metric_name="evenness", kingdom="bacteria")
input_params[["bacteria_phylo_d"]] <- list(alpha_div_in=bacteria_phylo_d, metric_name="phylo_d", kingdom="bacteria")
input_params[["bacteria_shannon"]] <- list(alpha_div_in=bacteria_shannon, metric_name="shannon", kingdom="bacteria")
input_params[["fungi_richness"]] <- list(alpha_div_in=fungi_richness, metric_name="richness", kingdom="fungi")
input_params[["fungi_evenness"]] <- list(alpha_div_in=fungi_evenness, metric_name="evenness", kingdom="fungi")
input_params[["fungi_phylo_d"]] <- list(alpha_div_in=fungi_phylo_d, metric_name="phylo_d", kingdom="fungi")
input_params[["fungi_shannon"]] <- list(alpha_div_in=fungi_shannon, metric_name="shannon", kingdom="fungi")

pheno_spearman_out <- lapply(input_params, function(x) { root_cor_alpha_and_pheno(root_phenos_in = root_phenos,
                                                                                  alpha_div_in = x$alpha_div_in,
                                                                                  metric_name = x$metric_name,
                                                                                  kingdom = x$kingdom)})


pheno_spearman <- do.call(rbind, pheno_spearman_out)
rownames(pheno_spearman) <- NULL

write.table(x = pheno_spearman, file = "alpha_div_cor/pheno_spearman_alpha_cor.tsv", col.names = TRUE, row.names=FALSE, sep="\t", quote=FALSE)
