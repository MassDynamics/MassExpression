#' This function orchestrates imputation, normalization and the binary limmastatistics accross all experimental comparisons
#' 
#' @param experimentDesign Output from constructSummarizedExperiment
#' @export runLimmaPipeline
#' @import data.table
#' @importFrom stringr str_c


runLimmaPipeline <- function(IntensityExperiment){
  
  
  intensityDF <- prepare_prot_int(IntensityExperiment)
  
  # Create Median Normalized Measurements in each Condition/Replicate
  intensityDF[Imputed == F, log2NIntNorm := log2NInt - median(log2NInt), by = list(Condition,Replicate)]
  
  # Imputation
  intensityDF <- imputeLFQ(intensityDF, 
                         id_type = "ProteinId", 
                         int_type = "log2NIntNorm",
                         f_imputeStDev = 0.3,
                         f_imputePosition= 1.8)
  
  
  # RunId will be unique to a row wherease replicate may not
  intensityDF[, RunId := str_c(Condition, Replicate, sep = ".")]
  
  #Run LIMMA
  resultsQuant <- limmaStatsFun(ID_type = "ProteinId",
                                   int_type = "log2NIntNorm",
                                   condition_col_name = "Condition",
                                   run_id_col_name = "RunId",
                                   rep_col_name = "Replicate",
                                   funDT = intensityDF)
  
  # Add DE analysis results to rowData
  rowData(IntensityExperiment) = merge(rowData(IntensityExperiment), 
                                       resultsQuant, by = "ProteinId", all.x = T)
  
  
  
  list(IntensityExperiment=IntensityExperiment,IntensityDF=intensityDF)
}
