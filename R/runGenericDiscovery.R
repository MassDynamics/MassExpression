#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @export

runGenericDiscovery <- function(experimentDesign, proteinIntensities){
  
  # Create Data Rep
  IntensityExperiment <- constructSummarizedExperiment(experimentDesign = experimentDesign, 
                                                       proteinIntensities = proteinIntensities)

  
  # Get Binary Statistic Comparisons and Long for Protein Intensity
  results <- runLimmaPipeline(IntensityExperiment)
  IntensityExperiment <- results[[1]]
  
  return(IntensityExperiment)
  
  # Currently object not used
  # intensities <- results[[2]]
}
