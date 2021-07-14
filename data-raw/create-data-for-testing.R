################################
# Output data to use for testing 
library(MassExpression)
design <- mq_lfq_data$design
intensities <- mq_lfq_data$intensities
parameters <- input_gen$parameters
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
     file = file.path("../tests/data/mq_lfq_output.RData"))



