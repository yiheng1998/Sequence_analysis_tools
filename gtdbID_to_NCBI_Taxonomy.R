# This code converts GTDB IDs to corresponding NCBI Taxonomy information.
# It supports batch operations and accepts valid inputs including 
# a single ID (e.g., "GB_GCA_000010645.1") or a vector of multiple IDs 
# (e.g., c("GB_GCA_000010645.1", "GB_GCA_000017645.1")).


library(rentrez)
library(xml2)
library(dplyr)
library(tidyr)
library(stringr)



batch_query_taxonomy <- function(acc_list) {
  result <- list()  
  for (acc in acc_list) {
    pure_acc <- str_extract(acc, "G[CF]A_.*")
    # Get UID
    search_res <- entrez_search(
      db = "assembly",
      term = pure_acc,  
      retmax = 1
    )
    if (search_res$count == 0) {
      warning(paste(acc, "No record"))
      result[[acc]] <- NA 
      next
    }
    uid <- search_res$ids
    
    # Get taxonomy information
    tryCatch({
      Sys.sleep(0.5)
      assembly_info <- entrez_summary(db = "assembly", id = uid)
      tax_id <- assembly_info$taxid
      tax_xml <- entrez_fetch(db = "taxonomy", id = tax_id, rettype = "xml")
      xml_doc <- read_xml(tax_xml)
      taxon_nodes <- xml_find_all(xml_doc, "//LineageEx/Taxon")  
      species_name <- xml_text(xml_find_first(xml_doc, "//Taxon[ScientificName]/ScientificName"))
      
      taxonomy_list <- list(
        Domain = NA,
        Phylum = NA,
        Class = NA,
        Order = NA,
        Family = NA,
        Genus = NA,
        Species = species_name
      )
      
      for (node in taxon_nodes) {
        rank <- xml_text(xml_find_first(node, "./Rank"))
        name <- xml_text(xml_find_first(node, "./ScientificName"))
        if (rank == "domain") taxonomy_list$Domain <- name
        else if (rank == "phylum") taxonomy_list$Phylum <- name
        else if (rank == "class") taxonomy_list$Class <- name
        else if (rank == "order") taxonomy_list$Order <- name
        else if (rank == "family") taxonomy_list$Family <- name
        else if (rank == "genus") taxonomy_list$Genus <- name
      }
      
      taxonomy_df <- as.data.frame(taxonomy_list, stringsAsFactors = FALSE)
      result[[acc]] <- taxonomy_df  
      cat(acc, "succeeded\n")
      
    }, error = function(e) {
      warning(paste(acc, "failed：", e$message))
      result[[acc]] <- NA
    })
  }
  
  result_df <- do.call(rbind, lapply(names(result), function(acc) {
    if (all(is.na(result[[acc]]))) {
      data.frame(
        acc = acc,
        Domain = NA, Phylum = NA, Class = NA, Order = NA, Family = NA, Genus = NA, Species = NA,
        stringsAsFactors = FALSE
      )
    } else {
      cbind(acc = acc, result[[acc]], stringsAsFactors = FALSE)
    }
  }))
  
  return(result_df) 
}
