#' This function orchestrates imputation, normalization and the binary limma statistics accross all experimental comparisons
#' 
#' @param longIntensityDT Output from `preProcess`
#' @param metadataExperiment Output from `preProcess`
#' @param comparisonType Type of model / pairwise comparison. One of, "all", "oneVSall", "custom", "spline". 
#' Default to "all"  
#' @param conditionsDict Condtion dictionary.
#' @param fitSeparateModels logical. TRUE to fit separate limma models for each pairwise comparisons 
#' (e.g. filtering and `lmFit` are run separately by comparison).
#' @param returnDecideTestColumn logical. If TRUE the row data of the `CompleteIntensityExperiment` will contain the output from 
#' `limma::decideTests`. If FALSE a single model is run for all contrasts.
#' @param conditionSeparator string. String used to separate up and down condition in output. 
#' 
#' @export 
#' 
#' @import data.table
#' @importFrom stringr str_c


fitModelTwoConditions <- function(longIntensityDT, 
                                 metadataExperiment,
                                 comparisonType = "all",
                                 orderCondition = NULL, 
                                 baselineInpuLevel = NULL, 
                                 customComparisonsCond = NULL, 
                                 conditionsDict,
                                 fitSeparateModels, 
                                 returnDecideTestColumn, 
                                 conditionSeparator
){
  
  # RunId will be unique to a row wheraes replicate may not
  #TODO why I cannot just use SampleName as the runId? 
  # ok with one condition
  longIntensityDT[, RunId := str_c(Condition, Replicate, sep = ".")]
  
  if(comparisonType %in% c("all", "oneVSall", "custom")){
    condLevels <- metadataExperiment$experimentType$condition1Levels
    condName <- metadataExperiment$experimentType$condition1Name
    print(paste0("Running pairwise comparison type", comparisonType ," with condition: ", condName))
    
    pairwiseComp <- createPairwiseComparisons(comparisonType = comparisonType, 
                                              condLevels = condLevels, 
                                              conditionsDict = conditionsDict,
                                              orderCondition = orderCondition,
                                              baselineInpuLevel = baselineInpuLevel,
                                              customComparisonsCond = customComparisonsCond)
    
    
    print("Starting DE with limma...")
    resultsQuant <- limmaPairwiseOneCondition(featureIdType = "ProteinId",
                                              intensityType = "log2NIntNorm",
                                              conditionColname = "Condition",
                                              runIdColname = "RunId",
                                              repColname = "Replicate",
                                              longIntensityDT = longIntensityDT,
                                              metadataExperiment = metadataExperiment,
                                              fitSeparateModels = fitSeparateModels,
                                              pairwiseComparisons = pairwiseComp,
                                              returnDecideTestColumn=returnDecideTestColumn, 
                                              conditionSeparator=conditionSeparator)
    
    stats <- resultsQuant[["resultsModel"]]
    conditionComparisonMapping <- resultsQuant[["conditionComparisonMapping"]]
    
    print("Decode conditions")
    conditionComparisonMapping <- MassExpression:::condition_name_decode_comparison_mapping(dt = conditionComparisonMapping, 
                                                                                            dict=conditionsDict, 
                                                                                            conditionSeparator=conditionSeparator)
    longIntensityDT <- MassExpression:::condition_name_decode_intensity_data(dt=longIntensityDT, 
                                                                             dict=conditionsDict)
    stats <- MassExpression:::condition_name_decode_limma_table(dt=stats, dict=conditionsDict)
    
    print("Limma analysis completed.")
    
    list(limmaStats=stats,
         decodedLongIntensityDT = longIntensityDT,
         conditionComparisonMapping=conditionComparisonMapping)
    
  } else {
    print("Spline or other models not yet implemented")
    return(NULL)
  }
}