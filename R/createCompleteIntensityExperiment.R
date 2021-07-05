#' Adds `intensityDF` information to `IntensityExperiment` with statistics about 
#' the number of replicates and number of imputed values in each condition of interest
#' @param IntensityExperiment output from `constructSummarizedExperiment`
#' @param limmaStats Statistic table from limma output
#' @param normalisationAppliedToAssay one of 'None' or 'Median'
#' @param intensityDF Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @param conditionComparisonMapping a dataframe matching different limma statistics comparison strings to up and down conditions
#' @export createCompleteIntensityExperiment
#' 
#' @importFrom SummarizedExperiment rowData

createCompleteIntensityExperiment <- function(IntensityExperiment,
                                              limmaStats, 
                                              normalisationAppliedToAssay,
                                              longIntensityDT, 
                                              conditionComparisonMapping){
  
  # replace Intensity with missing values to normalized log scale with imputed values
  wideIntensityDT <- longIntensityDT %>% tidyr::pivot_wider(id_cols = "ProteinId", 
                                                     names_from = "SampleName", 
                                                     values_from = "log2NIntNorm")
  
  # Mask of imputed data
  wideImputedDT <- longIntensityDT %>% tidyr::pivot_wider(id_cols = "ProteinId", 
                                                            names_from = "SampleName", 
                                                            values_from = "Imputed")
  
  intensityMatrixWide <- wideIntensityDT %>% dplyr::select(-ProteinId)
  intensityMatrixWide <- as.matrix(intensityMatrixWide)
  rownames(intensityMatrixWide) <- wideIntensityDT$ProteinId
  
  imputedMatrixWide <- wideImputedDT %>% dplyr::select(-ProteinId)
  imputedMatrixWide <- as.matrix(imputedMatrixWide)
  rownames(imputedMatrixWide) <- wideImputedDT$ProteinId

  metadataComplete <- metadata(IntensityExperiment)
  metadataComplete$NormalisationAppliedToAssay <- normalisationAppliedToAssay
  
  CompleteIntensityExperiment <- SummarizedExperiment(assays = SimpleList(intensities = intensityMatrixWide, 
                                                                          imputedMask = imputedMatrixWide),
                                                      rowData = merge(rowData(IntensityExperiment),
                                                                      limmaStats, by = "ProteinId", all.x = T), 
                                                      colData = colData(IntensityExperiment),
                                                      metadata = metadataComplete)
    
  # add imputed and replicate counts to the final object
  imputedCounts <- computeImputedCounts(longIntensityDT)
  rowData(CompleteIntensityExperiment) <- S4Vectors::merge(rowData(CompleteIntensityExperiment),
                                                           as.data.frame(imputedCounts), by = "ProteinId", all.x = T)

  
  replicateCounts <- computeReplicateCounts(longIntensityDT)
  rowData(CompleteIntensityExperiment) <- S4Vectors::merge(rowData(CompleteIntensityExperiment),
                                                           as.data.frame(replicateCounts), by = "ProteinId", all.x = T)
  
  metadata(CompleteIntensityExperiment)$conditionComparisonMapping <- conditionComparisonMapping
  
  return(CompleteIntensityExperiment)
}
