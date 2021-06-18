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
int_table <- diann_load("data-s3/dia-nn/test-data-rpackage/data/diann_report.tsv")
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,],
                                             group.header="Protein.Group",
                                             id.header = "Precursor.Id",
                                             quantity.header = "Precursor.Normalised")
protein_ids <-  rownames(protein.groups)
assay <- as_data_frame(protein.groups)
assay$protein_id <- protein_ids
colnames(assay)[2:ncol(assay)] <- colnames(protein.groups)
sample_data <- tibble(runs = colnames(protein.groups), groups = c(1,1,0))

metadata <- tibble(species = "", measure_type = "", label_method = "DIA", software = "DIA-NN")

save(assay,sample_data,file = "./data/example_diann_input.rda")


################
# MaxQuant-SILAC
################

pg <- read_delim("data-s3/maxquant/SILAC/paper-sample-dataset/proteinGroups.txt",  "\t",
                 escape_double = FALSE, trim_ws = TRUE, col_types = cols(Reverse = col_character()))
summary <- read_delim("data-s3/maxquant/SILAC/paper-sample-dataset/summary.txt",  "\t",
                 escape_double = FALSE, trim_ws = TRUE)
design <- read_delim("data-s3/maxquant/SILAC/paper-sample-dataset/experimentalDesign.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
sample_data <- design %>% filter(Experiment %in% c("Heart", "COL", "PAN")) %>%
  unite(runs, Fraction, Experiment, sep = "_",remove = FALSE)

metadata <- tibble(species = "Mouse", measure_type = "normalised ratios", label_method = "SILAC", software = "MaxQuant")

protein.id.col <- "Majority protein IDs"
intensity_columns <- c("Ratio H/L normalized Heart", "Ratio H/L normalized PAN", "Ratio H/L normalized COL")
pg_silac <- pg[,c(protein.id.col,intensity_columns)]
assay <- pg_silac %>% rename(protein_id = protein.id.col,
                                "Heart" = `Ratio H/L normalized Heart`,
                                "PAN"  = `Ratio H/L normalized PAN`,
                                "COL" = `Ratio H/L normalized COL`)

save(assay,sample_data,metadata, file = "./data/example_silac_input.rda")

