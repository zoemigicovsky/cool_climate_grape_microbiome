### Accessory functions to be called from other scripts.

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
