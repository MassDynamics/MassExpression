################################
# Output data to use for testing 
library(MassExpression)
design <- mq_lfq_data$design
intensities <- mq_lfq_data$intensities
parameters <- mq_lfq_data$parameters
normalisation_method <- parameters[parameters[,1] == "UseNormalisationMethod",2]
species <- parameters[parameters[,1] == "Species",2]
labellingMethod <- parameters[parameters[,1] == "LabellingMethod",2]


listIntensityExperiments <- runGenericDiscovery(experimentDesign = design, 
                                                proteinIntensities = intensities, 
                                                normalisationMethod = normalisation_method, 
                                                species = species, 
                                                labellingMethod = labellingMethod)

expectedCompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
expectedIntensityExperiment <- listIntensityExperiments$IntensityExperiment
output_folder <- "tests/data"

expectedcomparisonExperiments <- 
  listComparisonExperiments(expectedCompleteIntensityExperiment)

save(expectedIntensityExperiment, 
     expectedCompleteIntensityExperiment, 
     expectedcomparisonExperiments,
     file = file.path("tests/data/mq_lfq_output.RData"))


##########################
# HER2 on Maxquant worflow
##########################

# test data created in dev/AQ-3.0-banchmark-MassExpression

#' Transform MaxQuant proteinGroups to universal input
MaxQuantTranform <- function(proteinGroups, design, species, useNormalisationMethod, labellingMethod){
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


###
library(stringr)

### Raw MaxQuant output
mq_data_input <- "../../data-s3/data-formats/maxquant/LFQ/HER2/data/"
protein.intensities <- read.csv(file.path(mq_data_input, "proteinGroups.txt"), sep = "\t", stringsAsFactors = F)
design <- read.csv(file.path(mq_data_input, "experimentDesign_original.txt"), sep = "\t", stringsAsFactors = F)

input_gen <- MaxQuantTranform(proteinGroups=protein.intensities, design = design, species = "human", useNormalisationMethod = "None", labellingMethod = "LFQ")
data_int_me <- input_gen$intensities
colnames(data_int_me) <- gsub("LFQ.", "lfq.", colnames(data_int_me))


## Load HER2 output after running LFQprocessing

library(LFQProcessing)
# Where the output from LFQ processing is saved
output_folder <- "../../data-s3/data-formats/maxquant/LFQ/HER2/data/output/"
# her2 <- protein_quant_runner(upload_folder, output_folder)
out_data <- read.delim(file.path(output_folder, "proteinGroups_quant.txt")) %>%
  dplyr::rename(Majority.protein.IDs=majority.protein.ids,
                Fasta.headers=fasta.headers)
out_data_1prot <- assign_uniprot_ids_mq(out_data)
out_data_bench <- out_data_1prot %>% dplyr::select("ProteinId",
                                                   starts_with("logFC."), 
                                                   starts_with("adj.P.Val."), 
                                                   starts_with("P.Value.")) 
colnames(out_data_bench)[2:ncol(out_data_bench)] <- paste0("Disco_", colnames(out_data_bench)[2:ncol(out_data_bench)])

# Intensities
lfq_intensities <- out_data_1prot %>% dplyr::select(ProteinId, starts_with("lfq."))

## Run MassExpression discovery
intensities <- input_gen$intensities
design <- input_gen$design
parameters <- input_gen$parameters
normalisation_method <- parameters[parameters[,1] == "UseNormalisationMethod",2]
species <- parameters[parameters[,1] == "Species",2]
labellingMethod <- parameters[parameters[,1] == "LabellingMethod",2]

results <- runGenericDiscovery(experimentDesign = design, 
                               proteinIntensities = intensities, 
                               normalisationMethod = normalisation_method, 
                               species = species, 
                               labellingMethod = labellingMethod)

completeExperiment <- results$CompleteIntensityExperiment
inten <- results$IntensityExperiment
comparisonExperiments <- 
  listComparisonExperiments(completeExperiment)

## Compare logFC and Pvals between MassExpression and LFQprocessing

expected_complete_new <- rowData(completeExperiment)
expected_comparison_new_run <- rowData(comparisonExperiments[[1]])

# out_data_bench : HER2 processed with LFQ processing

sum(!(out_data_bench$ProteinId %in% expected_complete_new$ProteinId))
sum(!(expected_complete_new$ProteinId %in% out_data_bench$ProteinId))
sum(!(expected_comparison_new_run$ProteinId %in% out_data_bench$ProteinId))

compare_new_run <- out_data_bench %>% left_join(as_tibble(expected_comparison_new_run))

expected_diff_fc_maxquant <- compare_new_run$FC - compare_new_run$Disco_logFC.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3

expected_diff_pval_maxquant <- compare_new_run$P.Value - compare_new_run$Disco_P.Value.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3

data_bench_maxquant <- out_data_bench

ggplot(compare_new_run, aes(x=FC, 
                            y = `Disco_logFC.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3`)) + 
  geom_point() +
  theme_bw() +
  geom_smooth()+
  ggtitle("Proteins used in pairwise comparison")



compare_new_run$NImpute <- compare_new_run$NImputed..Parental_SKBR3 + compare_new_run$NImputed..AZD8931_resistant_SKBR3_AZDRc
ggplot(compare_new_run, aes(x=-log10(P.Value), 
                            y = -log10(`Disco_P.Value.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3`), colour=NImpute)) + 
  geom_point() +
  theme_bw() +
  geom_smooth()+
  scale_colour_viridis_c()+
  ggtitle("using proteins used in pairwise comparison")


ggplot(compare_new_run, aes(x = Disco_logFC.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3, y = -log10(`Disco_P.Value.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3`))) + geom_point()

ggplot(compare_new_run, aes(x = FC, y = -log10(`P.Value`))) + geom_point()



save(data_bench_maxquant, expected_diff_fc_maxquant, expected_diff_pval_maxquant, file = "tests/data/HER2_maxquant_workflow.RData")



# Test data fro writeReplicatesData

SampleDT <- data.table(ProteinId = c(rep("Prot1", 10), rep("Prot2", 10)),
                       GeneName = NA, Description = NA,
                       log2NInt = rep(1.5, 20),
                       Condition = c(rep(c("A", "B"), each=5), rep(c("A", "B"), each=5)), 
                       Replicate = c(1:5, 1:5, 1:5, 1:5), 
                       Imputed =  c(0,0,0,0,0,
                                    0,0,0,1,1,
                                    1,0,0,0,0,
                                    1,1,0,0,0))

writeReplicateData(SampleDT, "tests/data")
