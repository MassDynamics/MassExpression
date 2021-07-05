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
  longIntensityDT <- normaliseIntensity(longIntensityDT=longIntensityDT,
                                        NormalisationMethod=NormalisationMethod)
  
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
  #rowDataNoStats <- rowData(IntensityExperiment)
  # rowData(IntensityExperiment) = merge(rowData(IntensityExperiment), 
  #                                     stats, by = "ProteinId", all.x = T)
  
  # SummarizedExperiment which contains the complete essay with imputed data and statistics
  # about the number of proteins imputed in each condition 
  CompleteIntensityExperiment <- createCompleteIntensityExperiment(IntensityExperiment,
                                                                   limmaStats=stats,
                                                                   normalisationAppliedToAssay = NormalisationMethod,
                                                                   longIntensityDT = longIntensityDT,
                                                                   conditionComparisonMapping = conditionComparisonMapping)
  
  list(IntensityExperiment=IntensityExperiment,
       CompleteIntensityExperiment=CompleteIntensityExperiment)
}




#' This function normalises (based on `NormalisationMethod`) the log2 intensities stored in the `log2NInt` column initialised by `initialiseLongIntensityDT`. 
#' 
#' @param longIntensityDT Output from `initialiseLongIntensityDT`
#' @param NormalisationMethod logical. Whether to perform median normalisation by condition. 
#' @export normaliseIntensity
#' @return log2NIntNorm If `NormalisationMethod` is not `None` a new column is added containing the `log2NInt` normalised intensities 
#' @details Normalisation is only applied to non imputed intensities. `log2NIntNorm` is NA for imputed intensities. 

normaliseIntensity <- function(longIntensityDT, NormalisationMethod){
  normaliseDT <- longIntensityDT
  # Create Median Normalized Measurements in each Condition/Replicate
  if (NormalisationMethod == "Median"){
    normaliseDT[Imputed == F, log2NIntNorm := log2NInt - median(log2NInt), by = list(Condition,Replicate)]
  }else if(NormalisationMethod == "None") {
    normaliseDT[, log2NIntNorm := log2NInt]
  } else {
    stop(paste0("Normalisation method ",NormalisationMethod," not availble.  Available methods are: `Median` and `None`"))
  }
  return(normaliseDT)
}