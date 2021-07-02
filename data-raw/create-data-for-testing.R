################################
# Output data to use for testing 
library(MassExpression)
design <- mq_lfq_data$design
intensities <- mq_lfq_data$intensities

listIntensityExperiments <- runGenericDiscovery(experimentDesign = design, 
                                                proteinIntensities = intensities, 
                                                NormalisationMethod = "None")

expectedCompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
expectedIntensityExperiment <- listIntensityExperiments$IntensityExperiment
output_folder <- "tests/data"

expectedcomparisonExperiments <- 
  listComparisonExperiments(CompleteIntensityExperiment)

save(expectedIntensityExperiment, 
     expectedCompleteIntensityExperiment, 
     expectedcomparisonExperiments,
     file = file.path(output_folder, "mq_lfq_output.RData"))


