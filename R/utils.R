
################
### protein intensities
################

#' This function converts the assay data in a protein intensity to a long data.table object with 
#' intensities for each protein Id and intensity column
#' 
#' @param IntensityExperiment Output from constructSummarizedExperiment
#' @export SEToLongDT
#' @importFrom SummarizedExperiment rowData colData assay
#'
#' @import data.table
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join as_tibble

SEToLongDT <- function(IntensityExperiment){
  wide <- data.frame(assay(IntensityExperiment))
  wide$ProteinId <- rowData(IntensityExperiment)$ProteinId
  long <- wide %>% pivot_longer(cols = all_of(colnames(assay(IntensityExperiment))), 
                                names_to = "IntensityColumn", 
                                values_to = "Intensity")
  
  long <- long %>% left_join(as_tibble(colData(IntensityExperiment)))
  long <- long %>% left_join(as_tibble(rowData(IntensityExperiment)))
  as.data.table(long)
}



#' This function performs the log2 conversion and initialise the imputed column
#' @param IntensityExperiment Output from constructSummarizedExperiment
#' @export initialiseLongIntensityDT

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


#' Adds `intensityDF` information to `IntensityExperiment` with statistics about 
#' the number of replicates and number of imputed values in each condition of interest
#' @param IntensityExperiment output from `constructSummarizedExperiment`
#' @param intensityDF Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @export createCompleteIntensityExperiment
#' 
#' @importFrom SummarizedExperiment rowData

createCompleteIntensityExperiment <- function(IntensityExperiment, longIntensityDT){
  
  CompleteIntensityExperiment <- IntensityExperiment
  
  # replace Intensity with missing values to normalized log scale with imputed values
  wideIntensityDT <- longIntensityDT %>% pivot_wider(id_cols = "ProteinId", 
                                               names_from = "IntensityColumn", 
                                               values_from = "log2NIntNorm")
  
  intensityMatrixWide <- wideIntensityDT %>% dplyr::select(-ProteinId)
  intensityMatrixWide <- as.matrix(intensityMatrixWide)
  rownames(intensityMatrixWide) <- wideIntensityDT$ProteinId
  
  assay(CompleteIntensityExperiment) <- intensityMatrixWide

  # add imputed and replicate counts to the final object
  imputedCounts <- computeImputedCounts(longIntensityDT)
  rowData(CompleteIntensityExperiment) <- merge(rowData(CompleteIntensityExperiment),
                                                           imputedCounts, by = "ProteinId", all.x = T)
  
  replicateCounts <- computeReplicateCounts(longIntensityDT)
  rowData(CompleteIntensityExperiment) <- merge(rowData(CompleteIntensityExperiment),
                                                           replicateCounts, by = "ProteinId", all.x = T)
  
  return(CompleteIntensityExperiment)
}