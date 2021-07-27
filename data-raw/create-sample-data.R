library(devtools)
# install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)
library(here)
library(readr)
library(dplyr)
library(tidyr)

print(here())

#' Transform MaxQuant proteinGroups to universal input
MaxQuantTransform <- function(proteinGroups, design, species, useNormalisationMethod, labellingMethod){
  proteinGroups <- convert_protein_groups_to_universal(proteinGroups)
  cols <- colnames(proteinGroups)[grepl("LFQ.intensity.", colnames(proteinGroups))]
  runs <- gsub("LFQ.intensity.", "", cols)
  groups <- gsub("_","",str_extract(runs,"_(.*)_"))
  
  experiment.design <- dplyr::as_tibble(cbind(cols,runs)) %>%
    dplyr::rename(mqExperiment = runs) %>%
    left_join(design)
  design <- experiment.design %>%
    dplyr::select(cols, experiment, techRep, bioRep) %>%
    dplyr::rename(SampleName = cols, Condition = experiment)
  intensities <- proteinGroups[,c(cols,"ProteinId")]
  
  parameters <- data.frame(X1 = c("Species", "UseNormalisationMethod", "LabellingMethod"),
                           X2 = c(species,useNormalisationMethod,labellingMethod))
  
  
  list(design=design, intensities=intensities, parameters=parameters)
}



#' This function removes bad rows in a maxquant protein intensity object
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



###########
## fragpipe
###########

# Download MaxQuant processed data from: https://app.massdynamics.com/experiments/1e2b17ed-835d-40bf-be75-c6a15f962ed0#/results-files
# Find raw data at PRIDE entry PXD026401

protein.intensities = read.csv("../../data-s3/data-formats/fragpipe/LFQ/PXD026401/data/combined_protein.tsv", sep = "\t", stringsAsFactors = F)
cols <- colnames(protein.intensities)[grepl("Total.Intensity", colnames(protein.intensities))]
runs <- gsub(".Total.Intensity", "", cols)
groups <- gsub("_[0-9]*$","",runs)
experiment.design <-as_tibble(cbind(cols,groups))
colnames(experiment.design) <- c("SampleName", "Condition")

protein.id.column = "Protein.ID"
gene.id.column = "Gene.Names"
description.id.column = "Description"
protein.intensities = protein.intensities[,c(protein.id.column,
                                             gene.id.column,
                                             description.id.column,
                                             experiment.design$SampleName)]
protein.intensities = protein.intensities %>% 
  rename(
    ProteinId = protein.id.column,
    GeneId = gene.id.column,
    Description = description.id.column
  )


design_fragpipe <- experiment.design
assay_fragpipe <- protein.intensities 

parameters <- data.frame(X1 = c("Species", "UseNormalisationMethod", "LabellingMethod"),
                                X2 = c("Human","None","LFQ"))

fragpipe_data <- list(design = design_fragpipe, intensities = assay_fragpipe, parameters = parameters)

write_delim(assay_fragpipe, "../universal-input-explore/data/fragpipe-lfq/intensities.tsv", delim = "\t")
save(fragpipe_data, file = "./data/example_fragpipe.rda")

######################
#### Maxquant LFQ HER2
######################

# Download MaxQuant processed data from: https://app.massdynamics.com/experiments/c6fc6c60-fe65-47cb-bd6d-a021f0ed8720#/results-files

experiment_home <- "../../data-s3/data-formats/maxquant/LFQ/HER2/data/"
protein.intensities <- read.csv(file.path(experiment_home, "proteinGroups.txt"), sep = "\t", stringsAsFactors = F)
design <- read.csv(file.path(experiment_home, "experimentDesign_original.txt"), sep = "\t", stringsAsFactors = F)

# debugonce(MaxQuantTranform(protein.intensities, design))

mq_lfq_data <- MaxQuantTranform(protein.intensities, design, "Human", "None", "LFQ")

write_delim(protein.intensities, "../universal-input-explore/data/maxquant-lfq/intensities.tsv", delim = "\t")
save(mq_lfq_data, file = "./data/example_mq_lfq.rda")


######
### PD
######
experiment_home = "../../data-s3/data-formats/proteome-discoverer/LFQ/Spectral_Clustering_LFQ/data/data/pd_complete/final_pd_node/"

protein.intensities = read.csv(file.path(experiment_home, "iPRG-no_clustering-no_mbr_Proteins.txt"), sep = "\t", stringsAsFactors = F)

cols <- colnames(protein.intensities)[grepl("apQuant.Area.", colnames(protein.intensities))]
runs <- gsub("apQuant.Area.", "", cols)
groups <- str_extract(runs,"sample[0-9]")

experiment.design <- as_tibble(cbind(cols,groups))
colnames(experiment.design) <- c("SampleName", "Condition")
experiment.design

extract_genename <- function(description){
  desl <- strsplit(description, split=" ")[[1]]
  gn <- desl[sapply(desl, function(z) str_detect(z, "GN="))]
  gn <- str_remove(gn, "GN=")
  if(length(gn)==0){
    gn <- list("")
  }
  return(gn)
}

protein.id.column = "Accession"
description.id.column = "Description"

protein.intensities$GeneName <- do.call(c, sapply(protein.intensities$Description, function(x) extract_genename(x)))

protein.intensities = protein.intensities[,c(protein.id.column,
                                             description.id.column,
                                             "GeneName",
                                             experiment.design$SampleName)]

protein.intensities = protein.intensities %>% 
  rename(
    ProteinId = protein.id.column,
    Description = description.id.column
  )

parameters <- data.frame(X1 = c("Species", "UseNormalisationMethod", "LabellingMethod"),
                                X2 = c("Yeast","None","LFQ"))


write_delim(protein.intensities, "../data/sample-data-app/generic-importer-input/protein_intensities.tsv", delim = "\t")
write_delim(parameters, "../data/sample-data-app/generic-importer-input/parameters.tsv", delim = "\t", col_names=FALSE)
write_delim(experiment.design, "../data/sample-data-app/generic-importer-input/experimental_design.tsv", delim = "\t", col_names=TRUE)

pd_iPRG_data <- list(design = experiment.design, intensities=protein.intensities, parameters=parameters)

save(pd_iPRG_data, file = "./data/example_pd_iRPG.rda")


