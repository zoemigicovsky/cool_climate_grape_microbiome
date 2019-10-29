### Accessory functions to be called from other scripts.

library(VennDiagram)

qiime2_taxa_breakdown <- function(taxa_in) {
  # Add taxonomic levels to QIIME2 taxa table.
  
  taxa_in$species <- taxa_in$taxonomy
  
  # Remove "Amiguous_taxa" and replace "uncultured" labels with Unclassified.
  taxa_in$species <- gsub(";Ambiguous_taxa", "", taxa_in$species)
  taxa_in$species <- gsub("metagenome", "Unclassified", taxa_in$species)
  taxa_in$species <- gsub("uncultured[^;]*", "Unclassified", taxa_in$species)
  
  # Add in missing taxa labels as Unclassified.
  for(i in 0:6) {
    
    label = paste("D_", as.character(i), "__", sep="") 
    
    rows_missing_label <- grep(label, taxa_in$species, invert = TRUE)
    
    unclassified_label <- paste(label, "Unclassified", sep="")
    
    taxa_in$species[rows_missing_label] <- paste(taxa_in$species[rows_missing_label],
                                                 unclassified_label,
                                                 sep=";")
    
  }
  
  taxa_in$kingdom <- gsub(";D_1__.*$", "", taxa_in$species)
  taxa_in$phylum <- gsub(";D_2__.*$", "", taxa_in$species)
  taxa_in$class <- gsub(";D_3__.*$", "", taxa_in$species)
  taxa_in$order <- gsub(";D_4__.*$", "", taxa_in$species)
  taxa_in$family <- gsub(";D_5__.*$", "", taxa_in$species)
  taxa_in$genus <- gsub(";D_6__.*$", "", taxa_in$species)
  
  return(taxa_in)
  
}


UNITE_qiime2_taxa_breakdown <- function(taxa_in) {
  # Add taxonomic levels to QIIME2 taxa table.
  # This function expects a taxonomy table with the UNITE database ids.
  
  taxa_in$species <- taxa_in$taxonomy
  
  # Remove "Amiguous_taxa" and replace "uncultured" labels with Unclassified.
  taxa_in$species <- gsub(";Ambiguous_taxa", "", taxa_in$species)
  taxa_in$species <- gsub("Unassigned", "", taxa_in$species)
  taxa_in$species <- gsub("metagenome", "Unclassified", taxa_in$species)
  taxa_in$species <- gsub("uncultured[^;]*", "Unclassified", taxa_in$species)
  taxa_in$species <- gsub("unidentified[^;]*", "Unclassified", taxa_in$species)
  
  tax_labels <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Add in missing taxa labels as Unclassified.
  for(label in tax_labels) {
    
    rows_missing_label <- grep(label, taxa_in$species, invert = TRUE)
    
    unclassified_label <- paste(label, "Unclassified", sep="")
    
    if(label != "k__") {
      taxa_in$species[rows_missing_label] <- paste(taxa_in$species[rows_missing_label],
                                                   unclassified_label,
                                                   sep=";")
    } else {
      taxa_in$species[rows_missing_label] <- unclassified_label
    }
    
  }
  
  taxa_in$kingdom <- gsub(";p__.*$", "", taxa_in$species)
  taxa_in$phylum <- gsub(";c__.*$", "", taxa_in$species)
  taxa_in$class <- gsub(";o__.*$", "", taxa_in$species)
  taxa_in$order <- gsub(";f__.*$", "", taxa_in$species)
  taxa_in$family <- gsub(";g__.*$", "", taxa_in$species)
  taxa_in$genus <- gsub(";s__.*$", "", taxa_in$species)
  
  return(taxa_in)
  
}


threeWayVennPercentGenus <- function(metadata, meta_col, asv_abun, taxa_df, meta_cat,
                                     labels=c("cat1", "cat2", "cat3"), colours=c("#009E73", "#E69F00", "#56B4E9")) {
  
  # Get sets of ASVs and genera present in each grouping (root of cover crops, grape roots, and grape soil).
  category1_samples <- rownames(metadata)[which(metadata[, meta_col] == meta_cat[1])]
  category2_samples <- rownames(metadata)[which(metadata[, meta_col] == meta_cat[2])]
  category3_samples <- rownames(metadata)[which(metadata[, meta_col] == meta_cat[3])]
  
  category1_ASVs <- rownames(asv_abun)[which(rowSums(asv_abun[, category1_samples]) > 0)]
  category1_genera <- unique(taxa_df[category1_ASVs, "genus"])
  
  category2_ASVs <- rownames(asv_abun)[which(rowSums(asv_abun[, category2_samples]) > 0)]
  category2_genera <- unique(taxa_df[category2_ASVs, "genus"])
  
  category3_ASVs <- rownames(asv_abun)[which(rowSums(asv_abun[, category3_samples]) > 0)]
  category3_genera <- unique(taxa_df[category3_ASVs, "genus"])
  
  category1_genera_count <- length(category1_genera)
  category2_genera_count <- length(category2_genera)
  category3_genera_count <- length(category3_genera)
  category1_2_genera_count <- length(which(category1_genera %in% category2_genera))
  category2_3_genera_count <- length(which(category2_genera %in% category3_genera))
  category1_3_genera_count <- length(which(category1_genera %in% category3_genera))
  category1_2_3_genera_count_TMP <- category1_genera[which(category1_genera %in% category2_genera)]
  category1_2_3_genera_count <- length(category1_2_3_genera_count_TMP[which(category1_2_3_genera_count_TMP %in% category3_genera)])
  
  grid.newpage()
  
  venn_out <- draw.triple.venn(area1=category1_genera_count,
                               area2=category2_genera_count,
                               area3=category3_genera_count,
                               n12 = category1_2_genera_count,
                               n23 = category2_3_genera_count,
                               n13 = category1_3_genera_count,
                               n123 = category1_2_3_genera_count,
                               category = labels,
                               scaled=TRUE,
                               fill = colours,
                               print.mode="percent",
                               cex=rep(2, 7),
                               cat.cex=rep(2, 3))
  
  return(venn_out)
  
}
