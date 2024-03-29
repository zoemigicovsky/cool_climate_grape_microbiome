# Run PERMANOVA based on rootstock, rootdepth, and the interaction based on several distance metrics.

rm(list=ls(all=TRUE))

# Set seed for reproducibility.
set.seed(20191021)

library(vegan)

# Initialize empty list for PERMANOVA output.
permanova_output <- list()
permanova_output[["bacteria"]] <- list()
permanova_output[["fungi"]] <- list()


# Calculate PERMANOVA for bacteria:
bacteria_bray_curtis <- as.dist(read.table("data/diversity_files/bacteria/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv",
                                           header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

bacteria_unweighted_unifrac <- as.dist(read.table("data/diversity_files/bacteria/diversity_grape/unweighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                                  header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

bacteria_weighted_unifrac <- as.dist(read.table("data/diversity_files/bacteria/diversity_grape/weighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                                  header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

meta_vine <- read.table("data/metadata/root_depth_bacteria_metadata_grape_vine.txt",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#Reduce metadata down to samples in the distance matrix
bacteria_meta <- meta_vine[colnames(as.matrix(bacteria_bray_curtis)), ]

permanova_output[["bacteria"]][["bray_curtis"]] <- adonis(bacteria_bray_curtis ~ rootstock + root_depth + vine + rootstock * root_depth,
                                                          data = bacteria_meta, permutations=9999)

permanova_output[["bacteria"]][["unweighted_unifrac"]] <- adonis(bacteria_unweighted_unifrac ~ rootstock + root_depth + vine + rootstock * root_depth,
                                                                 data = bacteria_meta, permutations=9999)

permanova_output[["bacteria"]][["weighted_unifrac"]] <- adonis(bacteria_weighted_unifrac ~ rootstock + root_depth + vine + rootstock * root_depth,
                                                               data = bacteria_meta, permutations=9999)


# Calculate PERMANOVA for fungi:

fungi_bray_curtis <- as.dist(read.table("data/diversity_files/fungi/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv",
                                           header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

fungi_unweighted_unifrac <- as.dist(read.table("data/diversity_files/fungi/diversity_grape/unweighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                                  header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

fungi_weighted_unifrac <- as.dist(read.table("data/diversity_files/fungi/diversity_grape/weighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                                header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

meta_vine <- read.table("data/metadata/root_depth_bacteria_metadata_grape_vine.txt",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

#Reduce metadata down to samples in the distance matrix
fungi_meta <- meta_vine[colnames(as.matrix(fungi_bray_curtis)), ]

permanova_output[["fungi"]][["bray_curtis"]] <- adonis(fungi_bray_curtis ~ rootstock + root_depth + vine + rootstock * root_depth,
                                                          data = fungi_meta, permutations=9999)

permanova_output[["fungi"]][["unweighted_unifrac"]] <- adonis(fungi_unweighted_unifrac ~ rootstock + root_depth + vine + rootstock * root_depth,
                                                                 data = fungi_meta, permutations=9999)

permanova_output[["fungi"]][["weighted_unifrac"]] <- adonis(fungi_weighted_unifrac ~ rootstock + root_depth + vine + rootstock * root_depth,
                                                               data = fungi_meta, permutations=9999)

saveRDS(object = permanova_output, file = "data/intermediate_RDS/betadiv_permanova_out.rds")


# Re-run PERMANOVA with rootstock coded as grafted vs ungrafted


# Initialize empty list for PERMANOVA output.
permanova_output <- list()
permanova_output[["bacteria"]] <- list()
permanova_output[["fungi"]] <- list()


# Calculate PERMANOVA for bacteria:
bacteria_bray_curtis <- as.dist(read.table("data/diversity_files/bacteria/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv",
                                           header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

bacteria_unweighted_unifrac <- as.dist(read.table("data/diversity_files/bacteria/diversity_grape/unweighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                                  header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

bacteria_weighted_unifrac <- as.dist(read.table("data/diversity_files/bacteria/diversity_grape/weighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                                header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

meta_vine <- read.table("data/metadata/root_depth_bacteria_metadata_grape_vine.txt",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

meta_vine$grafting <-  ifelse(meta_vine$rootstock=="new_york_muscat", "own", "grafted")

#Reduce metadata down to samples in the distance matrix
bacteria_meta <- meta_vine[colnames(as.matrix(bacteria_bray_curtis)), ]

permanova_output[["bacteria"]][["bray_curtis"]] <- adonis(bacteria_bray_curtis ~ grafting + root_depth + vine + grafting * root_depth,
                                                          data = bacteria_meta, permutations=9999)

permanova_output[["bacteria"]][["unweighted_unifrac"]] <- adonis(bacteria_unweighted_unifrac ~ grafting + root_depth + vine + grafting * root_depth,
                                                                 data = bacteria_meta, permutations=9999)

permanova_output[["bacteria"]][["weighted_unifrac"]] <- adonis(bacteria_weighted_unifrac ~ grafting + root_depth + vine + grafting * root_depth,
                                                               data = bacteria_meta, permutations=9999)

# Calculate PERMANOVA for fungi:

fungi_bray_curtis <- as.dist(read.table("data/diversity_files/fungi/diversity_grape/bray_curtis_distance_matrix_exported/distance-matrix.tsv",
                                        header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

fungi_unweighted_unifrac <- as.dist(read.table("data/diversity_files/fungi/diversity_grape/unweighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                               header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

fungi_weighted_unifrac <- as.dist(read.table("data/diversity_files/fungi/diversity_grape/weighted_unifrac_distance_matrix_exported/distance-matrix.tsv",
                                             header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1))

meta_vine <- read.table("data/metadata/root_depth_bacteria_metadata_grape_vine.txt",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

meta_vine$grafting <-  ifelse(meta_vine$rootstock=="new_york_muscat", "own", "grafted")

#Reduce metadata down to samples in the distance matrix
fungi_meta <- meta_vine[colnames(as.matrix(fungi_bray_curtis)), ]

permanova_output[["fungi"]][["bray_curtis"]] <- adonis(fungi_bray_curtis ~ grafting + root_depth + vine + grafting * root_depth,
                                                       data = fungi_meta, permutations=9999)

permanova_output[["fungi"]][["unweighted_unifrac"]] <- adonis(fungi_unweighted_unifrac ~ grafting + root_depth + vine + grafting * root_depth,
                                                              data = fungi_meta, permutations=9999)

permanova_output[["fungi"]][["weighted_unifrac"]] <- adonis(fungi_weighted_unifrac ~ grafting + root_depth + vine + grafting * root_depth,
                                                            data = fungi_meta, permutations=9999)

saveRDS(object = permanova_output, file = "data/intermediate_RDS/betadiv_grafting_permanova_out.rds")

