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
  
  print("Starting DE with limma...")

  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  # Create valid condition names
  encodedCondition <- condition_name_encoder(dt = longIntensityDT)
  longIntensityDT <- encodedCondition$dt
  conditionsDict <- encodedCondition$conditionsDict

  # Create Median Normalized Measurements in each Condition/Replicate
  longIntensityDT <- normaliseIntensity(longIntensityDT=longIntensityDT,
                                        normalisationMethod=normalisationMethod)

  # Imputation
  longIntensityDT <- imputeLFQ(longIntensityDT, 
                         id_type = "ProteinId", 
                         int_type = "log2NIntNorm",
                         f_imputeStDev = 0.3,
                         f_imputePosition= 1.8)

  
  # RunId will be unique to a row wheraes replicate may not
  longIntensityDT[, RunId := str_c(Condition, Replicate, sep = ".")]
  
  # Run LIMMA
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




#' This function normalises (based on `normalisationMethod`) the log2 intensities stored in the `log2NInt` column initialised by `initialiseLongIntensityDT`. 
#' 
#' @param longIntensityDT Output from `initialiseLongIntensityDT`
#' @param normalisationMethod logical. Whether to perform median normalisation by condition. 
#' @export normaliseIntensity
#' @return log2NIntNorm If `normalisationMethod` is not `None` a new column is added containing the `log2NInt` normalised intensities 
#' @details Normalisation is only applied to non imputed intensities. `log2NIntNorm` is NA for imputed intensities. 

normaliseIntensity <- function(longIntensityDT, normalisationMethod){
  normaliseDT <- longIntensityDT
  # Create Median Normalized Measurements in each Condition/Replicate
  if (normalisationMethod == "Median"){
    normaliseDT[Imputed == F, log2NIntNorm := log2NInt - median(log2NInt), by = list(Condition,Replicate)]
  }else if(normalisationMethod == "None") {
    normaliseDT[, log2NIntNorm := log2NInt]
  } else {
    stop(paste0("Normalisation method ",normalisationMethod," not availble.  Available methods are: `Median` and `None`"))
  }
  return(normaliseDT)
}


##########################
# Internal function to encode and decode condition names to avoi problems with makeContrats
##########################

#' Create sanitized condition names to be valid for input to makeContrasts
#' @param dt data.table with columns `Condition` to be encoded
#' @keywords internal
#' @noRd
condition_name_encoder <- function(dt) {
  dt_copy = copy(dt)
  
  conditionsDict <- data.table(original=dt_copy[, unique(Condition)])
  set.seed(255)
  conditionsDict[, safe := stringi::stri_rand_strings(.N, 5, pattern = "[A-Za-z]")]
  
  dt_copy[, original := Condition]
  dt_copy <- merge(dt_copy, conditionsDict, by="original", all = T)
  dt_copy[, Condition := safe]
  dt_copy[, `:=`(original = NULL, safe = NULL)]
  
  list(dt=dt_copy, conditionsDict=conditionsDict)
}


#' Decode condition names for intensity data to return initial names provided in input
#' @param dt data.table with columns `Condition` to be decoded
#' @param dict mapping dictionary to use for condition name decoding
#' @keywords internal
#' @noRd
condition_name_decode_intensity_data <- function(dt, dict){
  dt_copy = copy(dt)
  dt_copy <- merge(dt_copy, dict, by.x = "Condition", by.y = "safe", all.x = T)
  dt_copy[, `:=`(Condition = NULL)]
  setnames(dt_copy, "original", "Condition")
  dt_copy[, RunId := str_c(Condition, Replicate, sep = ".")]
  return(dt_copy)
}


#' Decode to initial condition names in limma table
#' @param dt limma table
#' @param dict mapping dictionary to use for condition name decoding
#' @keywords internal
#' @noRd
condition_name_decode_limma_table <- function(dt, dict){
  dt_copy <- copy(dt)
  dt_columns <- colnames(dt_copy)
  for (i in 1:dict[, .N]) {
    original <- dict[i, original]
    safe <- dict[i, safe]
    
    dt_columns <- str_replace(dt_columns, pattern = safe, replacement = original)
  }
  colnames(dt_copy) <- dt_columns
  return(dt_copy)
}

#' Decode to initial condition names in comparison mapping 
#' @param dt data.tabl with mapping of pairwise comparison performed
#' @param dict mapping dictionary to use for condition name decoding
#' @keywords internal
#' @noRd
condition_name_decode_comparison_mapping <- function(dt, dict, conditionSeparator){
  dt_copy <- copy(dt)
  for (i in 1:dict[, .N]) {
    original <- dict[i, original]
    safe <- dict[i, safe]
    
    dt_copy[, up.condition := str_replace(up.condition, pattern = safe, replacement = original)]
    dt_copy[, down.condition := str_replace(down.condition, pattern = safe, replacement = original)]
  }
  dt_copy[, comparison.string := paste(up.condition, down.condition, sep=conditionSeparator)]
  return(dt_copy)
}
