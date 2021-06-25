#########
# DIA-NN
#########

library(devtools)
# install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)
library(here)

print(here())

# This is the final report that once would get from DIA-NN.
# Using R package diann one can parse that output further to obtain peptides and protein groups.
# For the moment we assume that the user gave us a matrix of proteins inten and sample infos as input
df <- diann_load("../data-s3/dia-nn/test-data-rpackage/data/diann_report.tsv")
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,],
                                             group.header="Protein.Group",
                                             id.header = "Precursor.Id",
                                             quantity.header = "Precursor.Normalised")

protein.groups <- cbind(protein.groups, protein.groups)

protein_ids <-  rownames(protein.groups)
assay_diann <- as_data_frame(protein.groups)
assay_diann$ProteinId <- protein_ids
assay_diann <- assay_diann %>% separate(ProteinId, into = c("ProteinId"), sep = ";")
colnames(assay_diann)[1:(ncol(assay_diann)-1)] <- paste0(colnames(protein.groups),1:ncol(protein.groups))


sample_data_diann <- tibble(runs = paste0(colnames(protein.groups),1:ncol(protein.groups)), groups = c(1,1,0,1,0,0))


metadata <- tibble(species = "", measure_type = "", label_method = "DIA", software = "DIA-NN")

sample_data_diann <- sample_data_diann %>% dplyr::rename(Condition = groups, 
                                                         IntensityColumn = runs) %>%
  mutate(Replicate = paste0("F",1:nrow(sample_data_diann)))

diann_data <- list(intensities = assay_diann, design = sample_data_diann)
save(diann_data, file = "./data/example_diann.rda")


################
# MaxQuant-SILAC
################

pg <- read_delim("../data-s3/maxquant/SILAC/paper-sample-dataset/proteinGroups.txt",  "\t",
                 escape_double = FALSE, trim_ws = TRUE, col_types = cols(Reverse = col_character()))
summary <- read_delim("../data-s3/maxquant/SILAC/paper-sample-dataset/summary.txt",  "\t",
                 escape_double = FALSE, trim_ws = TRUE)
design <- read_delim("../data-s3/maxquant/SILAC/paper-sample-dataset/experimentalDesign.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
sample_data_mq_silac <- design %>% filter(Experiment %in% c("Heart", "COL", "PAN")) %>%
  unite(runs, Fraction, Experiment, sep = "_",remove = FALSE)

metadata_mq_silac <- tibble(species = "Mouse", measure_type = "normalised ratios", label_method = "SILAC", software = "MaxQuant")

protein.id.col <- "Majority protein IDs"
intensity_columns <- c("Ratio H/L normalized Heart", "Ratio H/L normalized PAN", "Ratio H/L normalized COL")
pg_silac <- pg[,c(protein.id.col,intensity_columns)]
assay_mq_silac <- pg_silac %>% rename(protein_id = all_of(protein.id.col),
                                "Heart" = `Ratio H/L normalized Heart`,
                                "PAN"  = `Ratio H/L normalized PAN`,
                                "COL" = `Ratio H/L normalized COL`)

mq_silac_data <- list(intensities = assay_mq_silac, design = sample_data_mq_silac, metadata = metadata_mq_silac)

save(mq_silac_data, file = "./data/example_mq_silac.rda")



###########
## fragpipe
###########

protein.intensities = read.csv("../data-s3/frag-pipe/combined_protein.tsv", sep = "\t", stringsAsFactors = F)
cols <- colnames(protein.intensities)[grepl("Total.Intensity", colnames(protein.intensities))]
runs <- gsub(".Total.Intensity", "", cols)
groups <- gsub("_[0-9]*$","",runs)
experiment.design <-as_data_frame(cbind(cols,runs,groups))
colnames(experiment.design) <- c("IntensityColumn", "Replicate", "Condition")

protein.id.column = "Protein.ID"
gene.id.column = "Gene.Names"
description.id.column = "Description"
protein.intensities = protein.intensities[,c(protein.id.column,
                                             gene.id.column,
                                             description.id.column,
                                             experiment.design$IntensityColumn)]
protein.intensities = protein.intensities %>% 
  rename(
    ProteinId = protein.id.column,
    GeneId = gene.id.column,
    Description = description.id.column
  )


design_fragpipe <- experiment.design
assay_fragpipe <- protein.intensities 

fragpipe_data <- list(intensities = assay_fragpipe, design = design_fragpipe)

save(fragpipe_data, file = "./data/example_fragpipe.rda")


#### Maxquant LFQ HER2
experiment_home <- "../data-s3/data-formats/maxquant/LFQ/HER2/data/"
protein.intensities <- read.csv(file.path(experiment_home, "proteinGroups.txt"), sep = "\t", stringsAsFactors = F)
design <- read.csv(file.path(experiment_home, "experimentDesign_original.txt"), sep = "\t", stringsAsFactors = F)

# debugonce(MaxQuantTranform(protein.intensities, design))
mq_lfq_data <- MaxQuantTranform(protein.intensities, design)

save(mq_lfq_data, file = "./data/example_mq_lfq.rda")

