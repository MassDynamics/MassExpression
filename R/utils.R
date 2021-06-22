
################
### protein intensities
################

#' This function converts the assay data in a protein intensity to a long data.table object with 
#' intensities for each protein Id and intensity column
#' 
#' @param IntensityExperiment Output from constructSummarizedExperiment
#' @export makeLongIntensityDF
#' @importFrom SummarizedExperiment rowData colData assay
#' @import data.table

makeLongIntensityDF <- function(IntensityExperiment){
  wide <- as.data.table(assay(IntensityExperiment))
  colnames(wide) <- IntensityExperiment$IntensityColumn
  wide$ProteinId <- rowData(IntensityExperiment)$ProteinId
  long <- melt(wide, id.vars = c("ProteinId"), variable.name = "IntensityColumn", value.name = "Intensity")
  long <- merge(long, SummarizedExperiment::colData(IntensityExperiment), by =  "IntensityColumn")
  as.data.table(long)
}


#' This function performs the log2 conversion and writes the imputed column
#' @param IntensityExperiment Output from constructSummarizedExperiment
#' @export prepareIntensityDF
#' @import data.table

prepareIntensityDF <- function(IntensityExperiment){
  protInt <- makeLongIntensityDF(IntensityExperiment)
  stopifnot(dim(protInt)[1]>0)
  protInt <- as.data.table(protInt)
  protInt[, Imputed := 0L]
  protInt[Intensity == 0, Imputed := 1L]
  protInt[, log2NInt := 0.0]
  protInt[Intensity > 0 , log2NInt := log2(Intensity)]
  
  protInt 
}


#' Count the number of imputed features by condition and write to wide table
#' @param intensityDF Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).  
#' @export computeImputedCounts
#' @importFrom dplyr group_by 
#' @importFrom stringr str_c

computeImputedCounts <- function(intensityDF){
  imputedCounts <- intensityDF %>% group_by(ProteinId, Condition) %>%
    summarize(NImputed = sum(Imputed))
  
  imputedCounts <- imputedCounts %>% pivot_wider(id_cols = "ProteinId",
                                                 names_from = "Condition",
                                                 values_from = "NImputed")
  
  colnames(imputedCounts)[-1] <- str_c("NImputed: ", colnames(imputedCounts)[-1])
  return(imputedCounts)
}

#' Count the replicates in each condition and write to wide table
#' @param intensityDF Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @export computeReplicateCounts
#' @importFrom dplyr group_by 
#' @importFrom stringr str_c

computeReplicateCounts <- function(intensityDF){
  replicateCounts <- intensityDF %>% group_by(ProteinId, Condition) %>%
    summarize(NReplicates = n())
  
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
#' @importFrom SummarizedExperiment rowData
#' @importFrom tidyr pivot_wider

createCompleteIntensityExperiment <- function(IntensityExperiment, intensityDF){
  
  CompleteIntensityExperiment <- IntensityExperiment
  
  # replace Intensity with missing values to normalized log scale with imputed values
  intensityDFWide <- intensityDF %>% pivot_wider(id_cols = "ProteinId", 
                                               names_from = "IntensityColumn", 
                                               values_from = "log2NIntNorm")
  
  intensityMatrixWide <- intensityDFWide %>% dplyr::select(-ProteinId)
  intensityMatrixWide <- as.matrix(intensityMatrixWide)
  rownames(intensityMatrixWide) <- intensityDFWide$ProteinId
  
  assay(CompleteIntensityExperiment) <- intensityMatrixWide

  # add imputed and replicate counts to the final object
  imputedCounts <- computeImputedCounts(intensityDF)
  rowData(CompleteIntensityExperiment) <- merge(rowData(CompleteIntensityExperiment),
                                                           imputedCounts, by = "ProteinId", all.x = T)
  
  replicateCounts <- computeReplicateCounts(intensityDF)
  rowData(CompleteIntensityExperiment) <- merge(rowData(CompleteIntensityExperiment),
                                                           replicateCounts, by = "ProteinId", all.x = T)
  
  return(CompleteIntensityExperiment)
}