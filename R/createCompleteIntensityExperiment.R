#' Adds `longIntensityDT` information to `IntensityExperiment` with statistics about 
#' the number of replicates and number of imputed values in each condition of interest
#' @param IntensityExperiment output from `createSummarizedExperiment`
#' @param limmaStats Statistic table from limma output
#' @param normalisationAppliedToAssay one of 'None' or 'Median'
#' @param longIntensityDT Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @param conditionComparisonMapping a dataframe matching different limma statistics comparison strings to up and down conditions

#' @export createCompleteIntensityExperiment
#' 
#' @importFrom SummarizedExperiment rowData

createCompleteIntensityExperiment <- function(IntensityExperiment,
                                              limmaStats, 
                                              normalisationAppliedToAssay,
                                              longIntensityDT, 
                                              conditionComparisonMapping){
  
  #########
  # 1. Create intensity matrix and mask of imputed data
  #########
  
  # replace Intensity with missing values to normalized log scale with imputed values
  wideIntensityDT <- longIntensityDT %>% tidyr::pivot_wider(id_cols = "ProteinId", 
                                                     names_from = "SampleName", 
                                                     values_from = "log2NIntNorm") %>%
    mutate(ProteinId = forcats::fct_reorder(ProteinId, rownames(IntensityExperiment)))
  
  # Mask of imputed data
  wideImputedDT <- longIntensityDT %>% tidyr::pivot_wider(id_cols = "ProteinId", 
                                                            names_from = "SampleName", 
                                                            values_from = "Imputed")%>%
    mutate(ProteinId = forcats::fct_reorder(ProteinId, rownames(IntensityExperiment)))
  
  intensityMatrixWide <- wideIntensityDT %>% dplyr::select(-ProteinId)
  intensityMatrixWide <- as.matrix(intensityMatrixWide)
  rownames(intensityMatrixWide) <- wideIntensityDT$ProteinId
  
  imputedMatrixWide <- wideImputedDT %>% dplyr::select(-ProteinId)
  imputedMatrixWide <- as.matrix(imputedMatrixWide)
  rownames(imputedMatrixWide) <- wideImputedDT$ProteinId

  # Order coldata
  intensityMatrixWide <- intensityMatrixWide[, colnames(IntensityExperiment)]
  imputedMatrixWide <- imputedMatrixWide[,colnames(IntensityExperiment)]
  
  #########
  # 2. Create Create metadata
  #########
  metadataComplete <- metadata(IntensityExperiment)
  metadataComplete$NormalisationAppliedToAssay <- normalisationAppliedToAssay
  
  #########
  # 3. Join limma stats, imputed and replicate counts
  #########
  imputedCounts <- computeImputedCounts(longIntensityDT)
  replicateCounts <- computeReplicateCounts(longIntensityDT)
  
  rowDataComplete <- data.frame(rowData(IntensityExperiment)) %>% left_join(limmaStats)
  rowDataComplete <- rowDataComplete %>% left_join(imputedCounts) 
  rowDataComplete <- rowDataComplete %>% left_join(replicateCounts)

  
  #########
  # 4. Create final summaized experiment
  #########
  # Validate data for summarised experiments
  .validateInputCompleteExperiment(rowDataComplete, IntensityExperiment, intensityMatrixWide, imputedMatrixWide)
  
  CompleteIntensityExperiment <- SummarizedExperiment(assays = SimpleList(intensities = intensityMatrixWide, 
                                                                          imputedMask = imputedMatrixWide),
                                                      rowData = rowDataComplete,
                                                      colData = colData(IntensityExperiment),
                                                      metadata = metadataComplete)

  
  metadata(CompleteIntensityExperiment)$conditionComparisonMapping <- conditionComparisonMapping
  
  return(CompleteIntensityExperiment)
}



.validateInputCompleteExperiment <- function(rowDataComplete, IntensityExperiment, intensityMatrixWide, imputedMatrixWide){
  
  protRowInt <- sum(rowDataComplete$ProteinId != rownames(intensityMatrixWide)) == 0
  protRowImp <- sum(rowDataComplete$ProteinId != rownames(imputedMatrixWide)) == 0
  colInt <- sum(colData(IntensityExperiment)$SampleName != colnames(intensityMatrixWide)) == 0
  colImp <- sum(colData(IntensityExperiment)$SampleName != colnames(imputedMatrixWide)) == 0
  
  if(!protRowInt | !protRowImp){
    stop("Cannot create CompleteIntensityExperiment: ProteinId ordering of intensities and protein features is not the same.")
  }
  
  if(!colInt | !colImp){
    stop("Cannot create CompleteIntensityExperiment: colnames order of intensities and coldata is not the same.")
  }
  
}
