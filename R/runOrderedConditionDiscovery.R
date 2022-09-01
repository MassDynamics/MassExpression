#' This function orchestrates the MassExpression workflow (could be called by a  workflow step)
#' 
#' @param experimentDesign data.frame. Experiment design provided in input by the user. Required columsn are: `SampleName` and `Condition`.
#' @param proteinIntensities data.frame. Wide matrix of intensities. Rows are proteins and columns are SampleNames. Required column: `ProteinId`. 
#' @param comparisonType Type of model / pairwise comparison. One of, "all", "oneVSall", "custom", "spline". 
#' Default to "all"  
#' @param orderConditionsList list. Each entry of the list provides the ordering of the relative condition. 
#' E.g. list(Condition = c("2h", "4h", "6h)), in this case "2h" is considered the first level of the condition. 
#' @param baselineConditionList list. Each entry of the list provides the baseline of the relative condition.
#' list(Condition = c("2h")). This is used to perform the `oneVsall` comparison type. If `NULL`, either the first 
#' entry of `orderConditions` is used or an automatic detection is performed (imperfect solution). 
#' @param customComparisonsList list of data.frames. When `comparisonType` is `custom` this field is required 
#' to establish the custom levels to compare for each condition.
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

#' @details Options for `comparisonType` arguments. If `comparisonType = "all"`, `orderConditionsList`/`baselineConditionList` 
#' and `customComparisonsList` are not considered. If `comparisonType = "custom"` then `customComparisonsList` is required. 
#' If `comparisonType = "oneVSall"`, `orderConditionsList`/`baselineConditionList` can be provided to establish the order 
#' of the comparisons. If they are provided, one level of the condition is automatically selected based on automatic ordering
#' selection. 

#' @export 

#' @import log4r

# This workflow assumes that data have already been pre-processed/normalised

runOrderedConditionDiscovery <- function(experimentDesign, proteinIntensities, 
                                comparisonType = "all", 
                                orderConditionsList = NULL, #list
                                baselineConditionList = NULL, # list
                                customComparisonsList = NULL, # data.frame with left/right levels
                                fitSeparateModels = FALSE,
                                returnDecideTestColumn = FALSE, 
                                conditionSeparator = " - "){
  
  
  info(MassExpressionLogger(), "IMPORT DATA")
  # Create Data Rep
  IntensityExperiment <- importData(experimentDesign = experimentDesign,
                                    proteinIntensities = proteinIntensities)
  # now the data can be used for plotting and qc
  print(IntensityExperiment)
  
  metadataExperiment <- metadata(IntensityExperiment)
  
  info(MassExpressionLogger(), "PRE-PROCESS DATA")
  #TODO: return summarised experiment instead of long
  preProcessedData <- preProcess(IntensityExperiment=IntensityExperiment)
  longIntensityDT <- preProcessedData$longIntensityDT
  conditionsDict <- preProcessedData$conditionsDict
  info(MassExpressionLogger(), "PRE-PROCESS DONE!")
  
  info(MassExpressionLogger(), "Codition safe encoding dictionary")
  print(as.data.frame(conditionsDict))
  
  # could create output for plotting after norm
  ## TODO: normaliseIntensity using SampleName instead of Cond+ repl could cause issues with strange characters?
  
  if(metadataExperiment$experimentType$conditionOnly){
    condName <- metadataExperiment$experimentType$condition1Name  
    info(MassExpressionLogger(), 
         paste0("RUN DIFFERENTIAL EXPRESSION ANALYSIS FOR ONE CONDITION ONLY: ", condName))
    
    # pass information only 
    resultsLimma <- fitModelOneCondition(longIntensityDT=longIntensityDT,
                                         metadataExperiment = metadataExperiment,
                                         comparisonType = comparisonType,
                                         orderCondition = orderConditionsList[[condName]], 
                                         baselineInpuLevel = baselineConditionList[[condName]], 
                                         customComparisonsCond = customComparisonsList[[condName]],
                                         conditionsDict = conditionsDict[[condName]],
                                         fitSeparateModels = fitSeparateModels,
                                         returnDecideTestColumn = returnDecideTestColumn,
                                         conditionSeparator = conditionSeparator)
    
    info(MassExpressionLogger(), "CREATING OUTPUT SUMMARIZED EXPERIMENT")
    results <- createResults(IntensityExperiment,
                             limmaStats=resultsLimma$limmaStats,
                             longIntensityDT = resultsLimma$decodedLongIntensityDT,
                             conditionComparisonMapping = resultsLimma$conditionComparisonMapping) 
    
  } else {
    print("Limma model for two conditions not implemented yet.")
    
    results <- IntensityExperiment
  }
  
  print("Workflow completed.")
  
  return(results)
  
}
