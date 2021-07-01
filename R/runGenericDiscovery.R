#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @param NormalisationMethod Normalisation method. One of "none" or "median". 
#' @return List of two SummarisedExperiment objects: `IntensityExperiment` 
#' containing the raw intensities and  `CompleteIntensityExperiment` including 
#' imputed intensities. 
#' @export

runGenericDiscovery <- function(experimentDesign, proteinIntensities, 
                                NormalisationMethod="None"){
  
  # Create Data Rep
  IntensityExperiment <- constructSummarizedExperiment(experimentDesign = experimentDesign, 
                                                       proteinIntensities = proteinIntensities)

  # Get Binary Statistic Comparisons and complete experiment containinig imputed Protein Intensity
  results <- runLimmaPipeline(IntensityExperiment,
                              NormalisationMethod=NormalisationMethod)
  
  CompleteIntensityExperiment <- results$CompleteIntensityExperiment
  IntensityExperiment <- results$IntensityExperiment
  
  return(results)

}
