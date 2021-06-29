
#' Save RData with SummarizedExperiment object and ProteinViz json
#' 
#' @param IntensityExperiment XX.
#' @param CompleteIntensityExperiment XX
#' @param output_folder XX
#' @export

saveOutput <- function(IntensityExperiment, CompleteIntensityExperiment, output_folder){
  dir.create(file.path(output_folder), recursive = TRUE, showWarnings = FALSE)
  
  comparisonExperiments <- 
    listComparisonExperiments(CompleteIntensityExperiment)
  
  save(IntensityExperiment, 
       CompleteIntensityExperiment, 
       comparisonExperiments,
       file = file.path(output_folder, "DiscoveryQuant.RData"))
  
  writeProteinViz(CompleteIntensityExperiment, output_folder)
}
