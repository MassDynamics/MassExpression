#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @export

runGenericDiscovery <- function(experimentDesign, proteinIntensities){
  
  # Create Data Rep
  IntensityExperiment <- constructSummarizedExperiment(experimentDesign = experimentDesign, 
                                                       proteinIntensities = proteinIntensities)

  # Get Binary Statistic Comparisons and complete experiment containinig imputed Protein Intensity
  results <- runLimmaPipeline(IntensityExperiment)
  
  return(results)
}
