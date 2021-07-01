
#' Transform MaxQuant proteinGroups to universal input
#' @export MaxQuantTranform
MaxQuantTranform <- function(proteinGroups, design, species, useNormalisationMethod, labellingMethod){
  proteinGroups <- convert_protein_groups_to_universal(proteinGroups)
  cols <- colnames(proteinGroups)[grepl("LFQ.intensity.", colnames(proteinGroups))]
  runs <- gsub("LFQ.intensity.", "", cols)
  groups <- gsub("_","",str_extract(runs,"_(.*)_"))
  
  experiment.design <- dplyr::as_tibble(cbind(cols,runs)) %>%
    dplyr::rename(mqExperiment = runs) %>%
    left_join(design)
  design <- experiment.design %>%
    dplyr::select(cols, experiment, name, techRep, bioRep) %>%
    dplyr::rename(IntensityColumn = cols, Condition = experiment, Replicate = name)
  intensities <- proteinGroups[,c(cols,"ProteinId")]
  
  parameters <- data.frame(X1 = c("Species", "UseNormalisationMethod", "LabellingMethod"),
                                  X2 = c(species,useNormalisationMethod,labellingMethod))
  
  
  list(design=design, intensities=intensities, parameters=parameters)
}



#' This function removes bad rows in a maxquant protein intensity object
#' @export filter_columns_mq
filter_columns_mq <- function(proteinGroups){
  cols.to.filter = c("Reverse", "Potential.contaminant", "Only.identified.by.site", "Contaminant")
  for (col in intersect(colnames(proteinGroups),cols.to.filter)){
    proteinGroups = proteinGroups[proteinGroups[[col]] != "+",]
    proteinGroups[[col]] = NULL
    
  }
  
  return(proteinGroups) 
}

#' This function uses regex to get the first uniprot assession and the
#' corresponding gene and description and adds them to a protein intensity table 
#' coming from maxquant
#' @export assign_uniprot_ids_mq
#' @importFrom stringr str_extract_all str_extract str_c
assign_uniprot_ids_mq <- function(proteinGroups){
  
  uniprot_pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  major_ids <- str_extract_all(proteinGroups$Majority.protein.IDs, uniprot_pattern)
  num_major_ids <- table(unlist(lapply(major_ids, length)))
  selected_ids <- unlist(lapply(major_ids, function(x){x[1]}))
  
  gene_patterns = str_c(selected_ids, ".*GN=(\\S*) ")
  string_id_to_gene = str_extract(proteinGroups$Fasta.headers, gene_patterns) # get the string for the right assession
  string_just_gene = str_extract(string_id_to_gene, "GN=(\\S*) ")
  genes = gsub("[GN=| ]","", string_just_gene)
  
  string_id_to_gene = str_extract(proteinGroups$Fasta.headers, gene_patterns) # get the string for the right assession
  string_description = str_extract(string_id_to_gene, " (.*) ")
  descriptions = gsub("GN.*","", string_description )
  
  
  proteinGroups$ProteinId <- selected_ids
  proteinGroups$GeneId <- genes
  proteinGroups$Description <- descriptions
  
  return(proteinGroups)
}

#' This is a wrapper for the mq parser functions that can be used to go
#' from maxquant data to protein intensity tables we can work with
#' @export convert_protein_groups_to_universal
convert_protein_groups_to_universal <- function(proteinGroups){
  proteinGroups = filter_columns_mq(proteinGroups)
  proteinGroups = assign_uniprot_ids_mq(proteinGroups)
  return(proteinGroups)
}



