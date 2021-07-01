library(devtools)
# install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)
library(here)

print(here())


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

protein.intensities = read.csv("../../data-s3/data-formats/fragpipe/LFQ/PXD026401/data/combined_protein.tsv", sep = "\t", stringsAsFactors = F)
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

write_delim(assay_fragpipe, "../universal-input-explore/data/fragpipe-lfq/intensities.tsv", delim = "\t")
save(fragpipe_data, file = "./data/example_fragpipe.rda")


#### Maxquant LFQ HER2
experiment_home <- "../../data-s3/data-formats/maxquant/LFQ/HER2/data/"
protein.intensities <- read.csv(file.path(experiment_home, "proteinGroups.txt"), sep = "\t", stringsAsFactors = F)
design <- read.csv(file.path(experiment_home, "experimentDesign_original.txt"), sep = "\t", stringsAsFactors = F)

# debugonce(MaxQuantTranform(protein.intensities, design))
mq_lfq_data <- MaxQuantTranform(protein.intensities, design)

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

experiment.design <- as_tibble(cbind(cols,runs,groups))
colnames(experiment.design) <- c("IntensityColumn", "Replicate", "Condition")
experiment.design

protein.id.column = "Accession"
gene.id.column = NULL
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

write_delim(protein.intensities, "../universal-input-explore/data/pd-iRPG/intensities.tsv", delim = "\t")
pd_iPRG_data <- list(intensities=protein.intensities, design = experiment.design)

save(pd_iPRG_data, file = "./data/example_pd_iRPG.rda")


