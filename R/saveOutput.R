
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


#' Write limma statistics for users to read
#' @param CompleteIntensityExperiment produced by runLimmaPipeline
#' @param outputFolder A path to a folder to write the statistics to.
#' @export

writeLimmaStatisticsTable <- function(CompleteIntensityExperiment, outputFolder){
  write.table(rowData(CompleteIntensityExperiment), sep = "\t",
              file = file.path(outputFolder, "protein_limma_statistics.tsv"))
}