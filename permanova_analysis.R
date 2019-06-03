library(vegan)
library(tidyverse)
bacteria_meta <- read.delim("bacteria/root_depth_bacteria_metadata_grape_vine.txt", row.names=1)

bacteria_bray <- read.delim("bacteria/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv", header = T, row.names=1)

#Reduce meta data down to samples in the distance matrix
bacteria_meta <- bacteria_meta[rownames(bacteria_bray),]


permanova_bacteria <- adonis(bacteria_bray ~ rootstock+root_depth+rootstock*root_depth+vine,
       data = bacteria_meta, permutations=999, method = "bray")

permanova_bacteria

#Call:
# adonis(formula = bacteria_bray ~ rootstock + root_depth + rootstock *      root_depth, data = bacteria_meta, permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#                       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# rootstock             2    0.9435 0.47175  2.2464 0.12467  0.001 ***
# root_depth            2    0.6514 0.32572  1.5510 0.08608  0.037 *  
# rootstock:root_depth  4    0.7230 0.18075  0.8607 0.09553  0.833    
# Residuals            25    5.2500 0.21000         0.69372           
# Total                33    7.5679                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Minor but significant effect of rootstock (R2=0.12467, P=0.001) and root depth (R2=0.08608, P=0.018) on bacterial microbiome. The effect of rootstock is (slightly) stronger than depth, however, that's confounded with location partly as well since the rootstocks are in different rows. 

fungi_meta <- read.delim("fungi/root_depth_fungi_metadata_grape.tsv", row.names=1)

fungi_bray <- read.delim("fungi/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv", header = T, row.names=1)

#Reduce meta data down to samples in the distance matrix
fungi_meta <- fungi_meta[rownames(fungi_bray),]

permanova_fungi <- adonis(fungi_bray ~ rootstock+root_depth+rootstock*root_depth,
                             data = fungi_meta, permutations=999, method = "bray")

permanova_fungi

# #Call:
# adonis(formula = fungi_bray ~ rootstock + root_depth + rootstock *      root_depth, data = fungi_meta, permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#                       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# rootstock             2    0.7548 0.37739  2.3417 0.22260  0.009 **
# root_depth            2    0.4147 0.20733  1.2865 0.12229  0.197   
# rootstock:root_depth  3    0.2874 0.09581  0.5945 0.08477  0.959   
# Residuals            12    1.9339 0.16116         0.57034          
# Total                19    3.3907                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Rootstock is a significant fator (R2=0.22260, p=0.009) but root depth is not for fungal microbiome. 
