library(devtools)
# install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)
library(here)
library(readr)
library(dplyr)
library(tidyr)

print(here())


################
# MaxQuant-SILAC - NOT INCLUDED IN THE PACKGE
################

pg <- read_delim("../../data-s3/data-formats/maxquant/SILAC/paper-sample-dataset/data/proteinGroups.txt",  "\t",
                 escape_double = FALSE, trim_ws = TRUE, col_types = cols(Reverse = col_character()))
summary <- read_delim("../../data-s3/data-formats/maxquant/SILAC/paper-sample-dataset/data/summary.txt",  "\t",
                 escape_double = FALSE, trim_ws = TRUE)
design <- read_delim("../../data-s3/data-formats/maxquant/SILAC/paper-sample-dataset/data/experimentalDesign.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
sample_data_mq_silac <- design %>% filter(Experiment %in% c("Heart", "COL", "PAN")) %>%
  unite(runs, Fraction, Experiment, sep = "_",remove = FALSE)


protein.id.col <- "Majority protein IDs"
intensity_columns <- c("Ratio H/L normalized Heart", "Ratio H/L normalized PAN", "Ratio H/L normalized COL")
pg_silac <- pg[,c(protein.id.col,intensity_columns)]
assay_mq_silac <- pg_silac %>% rename(protein_id = all_of(protein.id.col),
                                "Heart" = `Ratio H/L normalized Heart`,
                                "PAN"  = `Ratio H/L normalized PAN`,
                                "COL" = `Ratio H/L normalized COL`)

parameters <- data.frame(X1 = c("Species", "UseNormalisationMethod", "LabellingMethod"),
                                X2 = c("Human","None","SILAC"))

mq_silac_data <- list(intensities = assay_mq_silac, design = sample_data_mq_silac, parameters = parameters)

save(mq_silac_data, file = "./data/example_mq_silac.rda")



###########
## fragpipe
###########

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

