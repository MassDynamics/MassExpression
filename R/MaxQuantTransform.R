#' This function removes bad rows in a maxquant protein intensity object
#' @export filter_columns_mq
filter_columns_mq <- function(protein.intensities){
  cols.to.filter = c("Reverse", "Potential.contaminant", "Only.identified.by.site", "Contaminant")
  for (col in intersect(colnames(protein.intensities),cols.to.filter)){
    print(col)
    protein.intensities = protein.intensities[protein.intensities[[col]] != "+",]
    protein.intensities[[col]] = NULL
    
  }
  
  return(protein.intensities) 
}

#' This function uses regex to get the first uniprot assession and the
#' corresponding gene and description and adds them to a protein intensity table 
#' coming from maxquant
#' @export assign_uniprot_ids_mq
#' @importFrom stringr str_extract_all str_extract str_c
assign_uniprot_ids_mq <- function(protein.intensities){
  
  uniprot_pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  major_ids <- str_extract_all(protein.intensities$Majority.protein.IDs, uniprot_pattern)
  num_major_ids <- table(unlist(lapply(major_ids, length)))
  selected_ids <- unlist(lapply(major_ids, function(x){x[1]}))
  
  gene_patterns = str_c(selected_ids, ".*GN=(\\S*) ")
  string_id_to_gene = str_extract(protein.intensities$Fasta.headers, gene_patterns) # get the string for the right assession
  string_just_gene = str_extract(string_id_to_gene, "GN=(\\S*) ")
  genes = gsub("[GN=| ]","", string_just_gene)
  
  string_id_to_gene = str_extract(protein.intensities$Fasta.headers, gene_patterns) # get the string for the right assession
  string_description = str_extract(string_id_to_gene, " (.*) ")
  descriptions = gsub("GN.*","", string_description )
  
  
  protein.intensities$ProteinId <- selected_ids
  protein.intensities$GeneId <- genes
  protein.intensities$Description <- descriptions
  
  return(protein.intensities)
}

#' This is a wrapper for the mq parser functions that can be used to go
#' from maxquant data to protein intensity tables we can work with
#' @export convert_protein_groups_to_universal
convert_protein_groups_to_universal <- function(proteinGroups){
  protein.intensities = filter_columns_mq(proteinGroups)
  protein.intensities = assign_uniprot_ids_mq(protein.intensities)
  return(protein.intensities)
}

