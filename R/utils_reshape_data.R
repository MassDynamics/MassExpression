#' This function converts the assay data in a protein intensity to a long data.table object with 
#' intensities for each protein Id and intensity column
#' 
#' @param IntensityExperiment Output from createSummarizedExperiment
#' @export SEToLongDT
#' @importFrom SummarizedExperiment rowData colData assay
#'
#' @import data.table
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join as_tibble

SEToLongDT <- function(IntensityExperiment, assayName = "raw"){
  wide <- as_tibble(assays(IntensityExperiment)[[assayName]])
  wide$ProteinId <- rowData(IntensityExperiment)$ProteinId
  long <- wide %>% pivot_longer(cols = all_of(colnames(assay(IntensityExperiment))), 
                                names_to = "SampleName", 
                                values_to = "Intensity")
  
  long <- long %>% left_join(as_tibble(colData(IntensityExperiment)))
  as.data.table(long)
}

#' This function performs the log2 conversion for intensities larger than 0 and initialise the imputed column.
#' @param IntensityExperiment Output from createSummarizedExperiment

#' @noRd
#' @keywords internal

initialiseLongIntensityDT <- function(IntensityExperiment){
  protInt <- SEToLongDT(IntensityExperiment)
  stopifnot(dim(protInt)[1]>0)
  protInt <- as.data.table(protInt)
  protInt[, Imputed := 0L]
  protInt[Intensity == 0, Imputed := 1L]
  protInt[, log2NInt := 0.0]
  protInt[Intensity > 0 , log2NInt := log2(Intensity)]
  as.data.table(protInt) 
}


#' Pivot long data table from long to wide format

#' @noRd
#' @keywords internal
#' 
pivotDTLongToWide <- function(longDT, idCol, colNamesFrom, fillValuesFrom){
  castingFormula <- as.formula(str_c(idCol, " ~ ", colNamesFrom))
  wideDT <-
    dcast.data.table(longDT, castingFormula, value.var = fillValuesFrom)
  wideDT <- wideDT[str_order(get(idCol), numeric = T)]
  featureNames <- wideDT[, get(idCol)]
  runNames <- colnames(wideDT[, 2:ncol(wideDT)])
  runNames <- str_sort(runNames, numeric = TRUE)
  
  wideMatrix <- as.matrix(wideDT[, runNames, with = FALSE])
  rownames(wideMatrix) <- featureNames
  return(wideMatrix)
}

#' Count the number of imputed features by condition and write to wide table
#' @param intensityDF Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).  
#' @export computeImputedCounts
#' 
#' @importFrom dplyr group_by summarize 
#' @importFrom stringr str_c
#' @importFrom tidyr pivot_wider

computeImputedCounts <- function(intensityDF){
  imputedCounts <- intensityDF %>% group_by(ProteinId, Condition) %>%
    dplyr::summarize(NImputed = sum(Imputed))
  
  imputedCounts <- imputedCounts %>% pivot_wider(id_cols = "ProteinId",
                                                 names_from = "Condition",
                                                 values_from = "NImputed")
  
  colnames(imputedCounts)[-1] <- str_c("NImputed: ", colnames(imputedCounts)[-1])
  return(imputedCounts)
}

#' Count the replicates in each condition and write to wide table
#' @param intensityDF Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @export computeReplicateCounts

computeReplicateCounts <- function(intensityDF){
  replicateCounts <- intensityDF %>% group_by(ProteinId, Condition) %>%
    dplyr::summarize(NReplicates = dplyr::n())
  
  replicateCounts <- replicateCounts %>% pivot_wider(id_cols = "ProteinId",
                                                 names_from = "Condition",
                                                 values_from = "NReplicates")
  colnames(replicateCounts)[-1] <- str_c("NReplicates: ", colnames(replicateCounts)[-1])
  return(replicateCounts)
}



