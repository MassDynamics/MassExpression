#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @return List of two SummarisedExperiment objects: `IntensityExperiment` 
#' containing the raw intensities and  `CompleteIntensityExperiment` including 
#' imputed intensities. 
#' @export

runGenericDiscovery <- function(experimentDesign, proteinIntensities, 
                                normalise=FALSE){
  
  # Create Data Rep
  IntensityExperiment <- constructSummarizedExperiment(experimentDesign = experimentDesign, 
                                                       proteinIntensities = proteinIntensities)

  # Get Binary Statistic Comparisons and complete experiment containinig imputed Protein Intensity
  results <- runLimmaPipeline(IntensityExperiment, normalise)
  
  CompleteIntensityExperiment <- results$CompleteIntensityExperiment
  IntensityExperiment <- results$IntensityExperiment
  
  return(results)

}
