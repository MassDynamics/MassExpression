#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @param normalisationMethod Normalisation method. One of "none" or "median". 
#' @param species Species. One of 'Human', 'Mouse', 'Yeast', 'Other'
#' @param labellingMethod One of 'LFQ' or 'TMT'
#' 
#' @return List of two SummarisedExperiment objects: `IntensityExperiment` 
#' containing the raw intensities and  `CompleteIntensityExperiment` including 
#' imputed intensities and the results of the limma DE analysis. 

#' @examples 
#' design <- fragpipe_data$design
#' intensities <- fragpipe_data$intensities
#' parameters <- fragpipe_data$parameters
#' normalisation_method <- parameters[parameters[,1] == "UseNormalisationMethod",2]
#' species <- parameters[parameters[,1] == "Species",2]
#' labellingMethod <- parameters[parameters[,1] == "LabellingMethod",2]
#' listIntensityExperiments <- runGenericDiscovery(experimentDesign = design, 
#' proteinIntensities = intensities, 
#' normalisationMethod = normalisation_method,
#' species = species, 
#' labellingMethod = labellingMethod)


#' @export

runGenericDiscovery <- function(experimentDesign, proteinIntensities, 
                                normalisationMethod="None", species, labellingMethod){
  
  listMetadata <- list(Species = species, 
                       LabellingMethod = labellingMethod, 
                       NormalisationAppliedToAssay = "None")
  
  # Create Data Rep
  IntensityExperiment <- constructSummarizedExperiment(experimentDesign = experimentDesign, 
                                                       proteinIntensities = proteinIntensities,
                                                       listMetadata = listMetadata)

  # Get Binary Statistic Comparisons and complete experiment containinig imputed Protein Intensity
  results <- runLimmaPipeline(IntensityExperiment,
                              normalisationMethod=normalisationMethod)
  
  CompleteIntensityExperiment <- results$CompleteIntensityExperiment
  IntensityExperiment <- results$IntensityExperiment
  
  return(results)

}
