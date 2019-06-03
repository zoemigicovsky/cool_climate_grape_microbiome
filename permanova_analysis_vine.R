library(vegan)
library(tidyverse)
bacteria_meta <- read.delim("bacteria/root_depth_bacteria_metadata_grape_vine.txt", row.names=1)

bacteria_bray <- read.delim("bacteria/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv", header = T, row.names=1)

#Reduce meta data down to samples in the distance matrix
bacteria_meta <- bacteria_meta[rownames(bacteria_bray),]


permanova_bacteria <- adonis(bacteria_bray ~ rootstock+root_depth+rootstock*root_depth+vine,
       data = bacteria_meta, permutations=999, method = "bray")

permanova_bacteria

# Call:
#   adonis(formula = bacteria_bray ~ rootstock + root_depth + rootstock *      root_depth + vine, data = bacteria_meta, permutations = 999,      method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# rootstock             2    0.9435 0.47175  2.7809 0.12467  0.001 ***
#   root_depth            2    0.6514 0.32572  1.9200 0.08608  0.006 ** 
#   vine                  9    2.5282 0.28091  1.6559 0.33407  0.001 ***
#   rootstock:root_depth  4    0.7305 0.18263  1.0765 0.09653  0.320    
# Residuals            16    2.7143 0.16964         0.35866           
# Total                33    7.5679                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#There are fewer samples in fungi so I can just use the same metadata table

fungi_meta <- bacteria_meta

fungi_bray <- read.delim("fungi/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv", header = T, row.names=1)

#Reduce meta data down to samples in the distance matrix
fungi_meta <- fungi_meta[rownames(fungi_bray),]

permanova_fungi <- adonis(fungi_bray ~ rootstock+root_depth+rootstock*root_depth+vine,
                             data = fungi_meta, permutations=999, method = "bray")

permanova_fungi
# 
# Call:
#   adonis(formula = fungi_bray ~ rootstock + root_depth + rootstock *      root_depth + vine, data = fungi_meta, permutations = 999,      method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#                       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# rootstock             2    0.7548 0.37739 2.59408 0.22260  0.006 **
#   root_depth            2    0.4147 0.20733 1.42518 0.12229  0.147   
# vine                  6    0.9488 0.15813 1.08697 0.27982  0.401   
# rootstock:root_depth  3    0.3996 0.13321 0.91567 0.11786  0.601   
# Residuals             6    0.8729 0.14548         0.25743          
# Total                19    3.3907                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Rootstock is a significant fator (R2=0.22260, p=0.009) but root depth is not for fungal microbiome. 
