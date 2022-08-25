#' This function orchestrates imputation, normalization and the binary limma statistics accross all experimental comparisons
#' 
#' @param IntensityExperiment Output from createSummarizedExperiment
#' @param normalisationMethod logical. Whether to perform median normalisation by condition. 
#' @param fitSeparateModels logical. TRUE to fit separate limma models for each pairwise comparisons 
#' (e.g. filtering and `lmFit` are run separately by comparison).
#' @param returnDecideTestColumn logical. If TRUE the row data of the `CompleteIntensityExperiment` will contain the output from 
#' `limma::decideTests`. If FALSE a single model is run for all contrasts.
#' @param conditionSeparator string. String used to separate up and down condition in output. 
#' 
#' @export runLimmaPipeline
#' @import data.table
#' @importFrom stringr str_c

runLimmaPipeline <- function(IntensityExperiment, 
                             normalisationMethod, 
                             fitSeparateModels, 
                             returnDecideTestColumn, 
                             conditionSeparator){
  
  print("Starting pre-processing")

  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  # Create valid condition names
  encodedCondition <- condition_name_encoder(dt = longIntensityDT, condNames = "Condition")
  longIntensityDT <- encodedCondition$dt
  conditionsDict <- encodedCondition$conditionsDict$Condition

  # Create Median Normalized Measurements in each Condition/Replicate
  longIntensityDT <- normaliseIntensity(longIntensityDT=longIntensityDT,
                                        normalisationMethod=normalisationMethod)

  # Imputation
  longIntensityDT <- imputeLFQ(longIntensityDT, 
                         id_type = "ProteinId", 
                         int_type = "log2NIntNorm",
                         f_imputeStDev = 0.3,
                         f_imputePosition = 1.8)

  
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

  print("Limma analysis completed. Creating output summarized experiment...")
  
  # SummarizedExperiment which contains the complete essay with imputed data and statistics
  # about the number of proteins imputed in each condition 
  CompleteIntensityExperiment <- createCompleteIntensityExperiment(IntensityExperiment,
                                                                   limmaStats=stats,
                                                                   normalisationAppliedToAssay = normalisationMethod,
                                                                   longIntensityDT = longIntensityDT,
                                                                   conditionComparisonMapping = conditionComparisonMapping)

  list(IntensityExperiment=IntensityExperiment,
       CompleteIntensityExperiment=CompleteIntensityExperiment, 
       longIntensityDT=longIntensityDT)
}





