#' @export

writeOutput <- function(IntensityExperiment, outputFolder){
  # Write Discovery Outputs
  write_protein_viz(outputFolder, IntensityExperiment)
}