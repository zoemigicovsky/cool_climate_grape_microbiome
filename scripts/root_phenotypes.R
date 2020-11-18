rm(list=ls(all.names=TRUE))

library(tidyverse)
library(broom)

setwd("/home/gavin/github_repos/root_depth/data/")

root_phenos <- read_tsv("root_phenotypes.txt") %>% filter(!is.na(sample_id))

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

