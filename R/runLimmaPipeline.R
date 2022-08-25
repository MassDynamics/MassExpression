#' This function orchestrates imputation, normalization and the binary limma statistics accross all experimental comparisons
#' 
#' @param longIntensityDT Output from `preProcess`
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

runLimmaPipeline <- function(longIntensityDT, 
                             conditionsDict,
                             fitSeparateModels, 
                             returnDecideTestColumn, 
                             conditionSeparator){

  # RunId will be unique to a row wheraes replicate may not
  #TODO why I cannot just use SampleName as the runId? 
  longIntensityDT[, RunId := str_c(Condition, Replicate, sep = ".")]
  
  # Run LIMMA
  print("Starting DE with limma...")
  resultsQuant <- limmaStatsFun(ID_type = "ProteinId",
                                   int_type = "log2NIntNorm",
                                   condition_col_name = "Condition",
                                   run_id_col_name = "RunId",
                                   rep_col_name = "Replicate",
                                   funDT = longIntensityDT,
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
}





