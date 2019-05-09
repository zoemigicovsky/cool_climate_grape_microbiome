library(tidyverse)
library(broom)
#Bacteria correlation with root phenotypes 

root_phenos <- read_tsv("root_phenotypes.txt") %>% filter(!is.na(sample_id))

#Alpha diversity

alpha_diversity <- read_tsv("bacteria/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_bacteria <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_bacteria <- root_phenos_bacteria %>% rename(trunk_diameter='trunk_diameter_at ground')

#Cannot compute exact p-values with ties, but they aren't significant anyway (even before multiple testing correction so it doesn't matter)
tidy(cor.test(root_phenos_bacteria$observed_otus,root_phenos_bacteria$trunk_diameter, method="spearman"))
#estimate statistic p.value method                          alternative
# -0.164     7619.   0.354 Spearman's rank correlation rho two.sided

tidy(cor.test(root_phenos_bacteria$observed_otus,root_phenos_bacteria$sec_root_mass, method="spearman"))
#estimate statistic p.value method                          alternative
#0.0587     6161.   0.742 Spearman's rank correlation rho two.sided  

tidy(cor.test(root_phenos_bacteria$observed_otus,root_phenos_bacteria$tert_root_mass, method="spearman"))
#estimate statistic p.value method                          alternative
#-0.0335     6764.   0.851 Spearman's rank correlation rho two.sided 

tidy(cor.test(root_phenos_bacteria$observed_otus,root_phenos_bacteria$sec_tert_root_mass, method="spearman")) 
#estimate statistic p.value method                          alternative
#0.0142     6452.   0.936 Spearman's rank correlation rho two.sided  

tidy(cor.test(root_phenos_bacteria$observed_otus,root_phenos_bacteria$total_root_mass, method="spearman"))
#estimate statistic p.value method                          alternative
# -0.0900     7134.   0.613 Spearman's rank correlation rho two.sided  

###now do the same thing with fungi

root_phenos <- read_tsv("root_phenotypes.txt") %>% filter(!is.na(sample_id))

#Alpha diversity

alpha_diversity <- read_tsv("fungi/diversity_grape/observed_otus_vector_exported/alpha-diversity.tsv") %>% rename(sample_id=X1)

#Reduce phenotype data down to samples with diversity metrics

root_phenos_fungi <- root_phenos %>% inner_join(alpha_diversity)
root_phenos_fungi <- root_phenos_fungi %>% rename(trunk_diameter='trunk_diameter_at ground')

#Cannot compute exact p-values with ties, but they aren't significant anyway (even before multiple testing correction so it doesn't matter)
tidy(cor.test(root_phenos_fungi$observed_otus,root_phenos_fungi$trunk_diameter, method="spearman"))
#estimate statistic p.value method                          alternative
# -0.264     1681.   0.261 Spearman's rank correlation rho two.sided  

tidy(cor.test(root_phenos_fungi$observed_otus,root_phenos_fungi$sec_root_mass, method="spearman"))
#estimate statistic p.value method                          alternative
# 0.0830     1220.   0.728 Spearman's rank correlation rho two.sided  

tidy(cor.test(root_phenos_fungi$observed_otus,root_phenos_fungi$tert_root_mass, method="spearman"))
#estimate statistic p.value method                          alternative
#0.144     1138.   0.544 Spearman's rank correlation rho two.sided  

tidy(cor.test(root_phenos_fungi$observed_otus,root_phenos_fungi$sec_tert_root_mass, method="spearman")) 
#estimate statistic p.value method                          alternative
#0.106     1188.   0.655 Spearman's rank correlation rho two.sided  

tidy(cor.test(root_phenos_fungi$observed_otus,root_phenos_fungi$total_root_mass, method="spearman"))
#estimate statistic p.value method                          alternative
# 0.0144     1311.   0.952 Spearman's rank correlation rho two.sided  