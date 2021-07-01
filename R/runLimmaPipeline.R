#' This function orchestrates imputation, normalization and the binary limma statistics accross all experimental comparisons
#' 
#' @param experimentDesign Output from constructSummarizedExperiment
#' @param NormalisationMethod logical. Whether to perform median normalisation by condition. 
#' @export runLimmaPipeline
#' @import data.table
#' @importFrom stringr str_c


runLimmaPipeline <- function(IntensityExperiment, NormalisationMethod){
  
  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  
  # Create Median Normalized Measurements in each Condition/Replicate
  if (NormalisationMethod == "Median"){
    longIntensityDT[Imputed == F, log2NIntNorm := log2NInt - median(log2NInt), by = list(Condition,Replicate)]
  }else if(NormalisationMethod == "None") {
    longIntensityDT[, log2NIntNorm := log2NInt]
  } else {
    stop(paste0("Normalisation method ",NormalisationMethod," not availble."))
  }
  
  # Imputation
  longIntensityDT <- imputeLFQ(longIntensityDT, 
                         id_type = "ProteinId", 
                         int_type = "log2NIntNorm",
                         f_imputeStDev = 0.3,
                         f_imputePosition= 1.8)
  
  
  # RunId will be unique to a row wherease replicate may not
  longIntensityDT[, RunId := str_c(Condition, Replicate, sep = ".")]
  
  # Run LIMMA
  resultsQuant <- limmaStatsFun(ID_type = "ProteinId",
                                   int_type = "log2NIntNorm",
                                   condition_col_name = "Condition",
                                   run_id_col_name = "RunId",
                                   rep_col_name = "Replicate",
                                   funDT = longIntensityDT)
  stats = resultsQuant[["stats"]]
  conditionComparisonMapping = resultsQuant[["conditionComparisonMapping"]]
  
  # data used for DE analysis - filtered, imputed, normalised
  # eset not used for now but would be good to produce stats QC related to limma: variance plots
  # eset = resultsQuant[["eset"]]
  
  # Add DE analysis results to rowData
  rowData(IntensityExperiment) = merge(rowData(IntensityExperiment), 
                                       stats, by = "ProteinId", all.x = T)
  
  # SummarizedExperiment which contains the complete essay with imputed data and statistics
  # about the number of proteins imputed in each condition 
  CompleteIntensityExperiment <- createCompleteIntensityExperiment(IntensityExperiment=IntensityExperiment,
                                                                   longIntensityDT = longIntensityDT,
                                                                   conditionComparisonMapping = conditionComparisonMapping)
  
  list(IntensityExperiment=IntensityExperiment,
       CompleteIntensityExperiment=CompleteIntensityExperiment)
}
