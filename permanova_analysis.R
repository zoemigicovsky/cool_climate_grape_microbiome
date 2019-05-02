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
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# rootstock   1   0.13270 0.13270 0.96585 0.08081  0.369  
# root_depth  2   0.54766 0.27383 1.99308 0.33351  0.020 *
#   Residuals   7   0.96174 0.13739         0.58568         
# Total      10   1.64210                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Root depth is a significant fator (R2=0.33351, p=0.020) but rootstock is not for fungal microbiome. 
