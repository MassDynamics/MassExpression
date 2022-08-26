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
#' @export runLimmaPipeline
#' @import data.table
#' @importFrom stringr str_c

# one condition

runLimmaPipeline <- function(longIntensityDT, 
                             metadataExperiment,
                             comparisonType = "all",
                             orderConditions = NULL, 
                             baselineCondition = NULL, 
                             customComparisons = NULL, 
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
    condName <- metadataExperiment$experimentType$condition1Name
    condLevels <- metadataExperiment$experimentType$condition1Levels
    print(paste0("Create pairwise comparisons for condition: ",condName))
  
    pairwiseComp <- createPairwiseComparisons(comparisonType = comparisonType, 
                                              condLevels = condLevels, 
                                              conditionsDict = conditionsDict,
                                              orderCondition = orderConditions[[condName]],
                                              baselineInpuLevel = baselineCondition[[condName]],
                                              customComparisonsCond = customComparisons[[condName]])

    
    print("Starting DE with limma...")
    
  # featureIdType = "ProteinId"
  # intensityType = "log2NIntNorm"
  # conditionColname = "Condition"
  # runIdColname = "RunId"
  # repColname = "Replicate"
  # longIntensityDT = longIntensityDT
  # pairwiseComparisons = pairwiseComp
  # returnDecideTestColumn=returnDecideTestColumn
  # conditionSeparator=conditionSeparator
  #   
    resultsQuant <- limmaPairwiseComparisons(featureIdType = "ProteinId",
                                     intensityType = "log2NIntNorm",
                                     conditionColname = "Condition",
                                     runIdColname = "RunId",
                                     repColname = "Replicate",
                                     longIntensityDT = longIntensityDT,
                                     metadataExperiment = metadataExperiment,
                                     pairwiseComparisons = pairwiseComp,
                                  returnDecideTestColumn=returnDecideTestColumn, 
                                  conditionSeparator=conditionSeparator)
    
    if(fitSeparateModels){
      stats <- resultsQuant[["statsSepModels"]]
    }else{
      stats <- resultsQuant[['statsOneModel']] 
    }
   
    conditionComparisonMapping <- resultsQuant[["conditionComparisonMapping"]]
    
    conditionComparisonMapping <- condition_name_decode_comparison_mapping(dt = conditionComparisonMapping, 
                                                                           dict=conditionsDict, 
                                                                           conditionSeparator=conditionSeparator)
    longIntensityDT <- condition_name_decode_intensity_data(dt=longIntensityDT, dict=conditionsDict)
    stats <- condition_name_decode_limma_table(dt=stats, dict=conditionsDict)
  
    print("Limma analysis completed.")
    
    list(limmaStats=stats,
         decodedLongIntensityDT = longIntensityDT,
         conditionComparisonMapping=conditionComparisonMapping)
    
  } else {
    print("Spline or other models not yet implemented")
    return(NULL)
  }
}