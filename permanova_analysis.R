library(vegan)
library(tidyverse)
bacteria_meta <- read.delim("bacteria/root_depth_bacteria_metadata_grape.tsv", row.names=1)

bacteria_bray <- read.delim("bacteria/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv", header = T, row.names=1)

#Reduce meta data down to samples in the distance matrix
bacteria_meta <- bacteria_meta[rownames(bacteria_bray),]


permanova_bacteria <- adonis(bacteria_bray ~ rootstock+root_depth,
       data = bacteria_meta, permutations=999, method = "bray")

permanova_bacteria

# Call:
#   adonis(formula = t(bacteria_bray) ~ rootstock + root_depth, data = bacteria_meta,      permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#             Df SumsOfSqs MeanSqs F.Model    R2  Pr(>F)    
# rootstock   2    0.9233 0.46166  2.4465 0.13145  0.001 ***
# root_depth  2    0.6286 0.31428  1.6654 0.08948  0.018 *  
# Residuals  29    5.4725 0.18871         0.77907           
# Total      33    7.0244                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Minor but significant effect of rootstock (R2=0.13145, P=0.001) and root depth (R2=0.08948, P=0.018) on bacterial microbiome. The effect of rootstock is (slightly) stronger than depth, however, that's confounded with location partly as well since the rootstocks are in different rows. 

fungi_meta <- read.delim("fungi/root_depth_fungi_metadata_grape.tsv", row.names=1)

fungi_bray <- read.delim("fungi/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv", header = T, row.names=1)

#Reduce meta data down to samples in the distance matrix
fungi_meta <- fungi_meta[rownames(fungi_bray),]

permanova_fungi <- adonis(fungi_bray ~ rootstock+root_depth,
                             data = fungi_meta, permutations=999, method = "bray")

permanova_fungi

#Call:
# adonis(formula = fungi_bray ~ rootstock + root_depth, data = fungi_meta,      permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# rootstock   2    0.7548 0.37739  2.5484 0.22260  0.002 **
# root_depth  2    0.4147 0.20733  1.4001 0.12229  0.119   
# Residuals  15    2.2213 0.14809         0.65511          
# Total      19    3.3907                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Rootstock is a significant fator (R2=0.22260, p=0.002) but root depth is not for fungal microbiome. 
