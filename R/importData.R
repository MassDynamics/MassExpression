#' This function creates a summarized experiment object from a protein intensity table and experiment design. 
#' 
#' @param experimentDesign data.frame. Experiment design provided in input by the user. Required columsn are: `SampleName` and `Condition`.
#' @param proteinIntensities data.frame. Wide matrix of intensities. Rows are proteins and columns are SampleNames. Required column: `ProteinId`. 
#' @param normalisationMethod Normalisation method. One of "None" or "Median". 
#' @param species Species. One of 'Human', 'Mouse', 'Yeast', 'Other'
#' @param labellingMethod One of 'LFQ' or 'TMT'
#' 
#' @return A SummarizedExperiment object
#' 
#' @export
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData

importData <- function(experimentDesign, proteinIntensities, 
                       normalisationMethod="None", species, 
                       labellingMethod){
  
  print("Sanitize import")
  # Sanitize invisible and Unicode characters for correct export
  experimentDesign <- sanitize_strings_in_dataframe(experimentDesign)
  proteinIntensities <- sanitize_strings_in_dataframe(proteinIntensities)
  
  print("Prepare list metadata")
  listMetadata <- list(Species = species,
                       LabellingMethod = labellingMethod, 
                       NormalisationAppliedToAssay = normalisationMethod)
  
  # Create Data Rep
  IntensityExperiment <- createSummarizedExperiment(experimentDesign = experimentDesign, 
                                                    proteinIntensities = proteinIntensities,
                                                    listMetadata = listMetadata)
  
  return(IntensityExperiment)
}
