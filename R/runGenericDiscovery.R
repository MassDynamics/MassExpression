#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign data.frame. Experiment design provided in input by the user. Required columsn are: `SampleName` and `Condition`.
#' @param proteinIntensities data.frame. Wide matrix of intensities. Rows are proteins and columns are SampleNames. Required column: `ProteinId`. 
#' @param normalisationMethod Normalisation method. One of "None" or "Median". 
#' @param species Species. One of 'Human', 'Mouse', 'Yeast', 'Other'
#' @param labellingMethod One of 'LFQ' or 'TMT'
#' @param fitSeparateModels logical. TRUE to fit separate limma models for each pairwise comparisons 
#' (e.g. filtering and `lmFit` are run separately by comparison).
#' @param returnDecideTestColumn logical. If TRUE the row data of the `CompleteIntensityExperiment` will contain the output from 
#' `limma::decideTests`. If FALSE a single model is run for all contrasts.
#' @param conditionSeparator string. String used to separate up and down condition in output. 
#' 
#' @return List of two SummarisedExperiment objects: `IntensityExperiment` 
#' containing the raw intensities and  `CompleteIntensityExperiment` including 
#' imputed intensities and the results of the limma DE analysis. 
#' 

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
                                normalisationMethod="None", species, 
                                labellingMethod, 
                                fitSeparateModels = TRUE,
                                returnDecideTestColumn = FALSE, 
                                conditionSeparator = " - "){
  
  
  print("IMPORT DATA")
  # Create Data Rep
  IntensityExperiment <- importData(experimentDesign = experimentDesign,
                                    proteinIntensities = proteinIntensities,
                                    normalisationMethod=normalisationMethod, 
                                    species = species, 
                                    labellingMethod)
  
  metadataExperiment <- metadata(IntensityExperiment)
  if(metadataExperiment$experimentType$conditionOnly){
    # Get Binary Statistic Comparisons and complete experiment containinig imputed Protein Intensity
    results <- runLimmaPipeline(IntensityExperiment,
                                normalisationMethod=normalisationMethod, 
                                fitSeparateModels=fitSeparateModels,
                                returnDecideTestColumn=returnDecideTestColumn, 
                                conditionSeparator = conditionSeparator)
    
  } else {
    print("Limma model for two conditions not implemented yet.")
    results <- IntensityExperiment
  }

  print("Workflow completed.")

  return(results)

}
