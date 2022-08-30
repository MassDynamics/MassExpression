#' Adds `longIntensityDT` information to `IntensityExperiment` with statistics about 
#' the number of replicates and number of imputed values in each condition of interest
#' @param IntensityExperiment output from `createSummarizedExperiment`
#' @param limmaStats Statistic table from limma output
#' @param longIntensityDT Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @param conditionComparisonMapping a dataframe matching different limma statistics comparison strings to up and down conditions

#' @export createCompleteIntensityExperiment
#' 
#' @importFrom SummarizedExperiment rowData

createCompleteIntensityExperiment <- function(IntensityExperiment,
                                              limmaStats,
                                              longIntensityDT, 
                                              conditionComparisonMapping){
  
  #########
  # 1. Create intensity matrix and mask of imputed data
  #########
  
  # # From wide to long matrix of counts
  intensityMatrixWide <- MassExpression:::pivotDTLongToWide(longIntensityDT,
                                                          idCol = "ProteinId",
                                                          colNamesFrom = "SampleName",
                                                          fillValuesFrom = "log2NIntNorm")
  intensityMatrixWide <- intensityMatrixWide[match(rownames(IntensityExperiment), rownames(intensityMatrixWide)),]

  ### Imputed value mask matrix ------
  imputedMatrixWide <- MassExpression:::pivotDTLongToWide(longIntensityDT,
                                                      idCol = "ProteinId",
                                                      colNamesFrom = "SampleName",
                                                      fillValuesFrom = "Imputed")
  imputedMatrixWide <- imputedMatrixWide[match(rownames(IntensityExperiment), rownames(imputedMatrixWide)),]
  
  # Order coldata
  intensityMatrixWide <- intensityMatrixWide[, colnames(IntensityExperiment)]
  imputedMatrixWide <- imputedMatrixWide[, colnames(IntensityExperiment)]
  
  #########
  # 2. Create Create metadata
  #########
  metadataComplete <- metadata(IntensityExperiment)
  
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
  MassExpression:::.validateInputCompleteExperiment(rowDataComplete, IntensityExperiment, intensityMatrixWide, imputedMatrixWide)
  
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
