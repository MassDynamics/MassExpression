library(readr)
input_folder <- "~/Desktop/input-universal/"
output_folder <- "~/Desktop/input-universal/"

design <- read_delim(file.path(input_folder, "design.tsv"), delim = "\t")
intensities <- read_delim(file.path(input_folder, "intensities.tsv"), delim = "\t")
# metadata <- read_delim(file.path(input_folder, "metadata.tsv"), delim = "\t")

# the flag should come from the metadata 
normalise_flag <- FALSE

# Workflow runner
runGenericDiscovery(experimentDesign = design, 
                    proteinIntensities = intensities,
                    normalise = normalise_flag, 
                    output_folder = output_folder,
                    save=FALSE)
