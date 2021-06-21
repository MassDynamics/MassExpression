
################
### protein intensities
################

#' This function converts the assay data in a protein intensity to a long data.table object with 
#' intensities for each protein Id and intensity column
#' 
#' @param IntensityExperiment Output from constructSummarizedExperiment
#' @export get_long_protein_intensity
#' @importFrom SummarizedExperiment rowData assay
#' @import data.table

get_long_protein_intensity <- function(IntensityExperiment){
  wide <- as.data.table(assay(IntensityExperiment))
  colnames(wide) <- IntensityExperiment$IntensityColumn
  wide$ProteinId <- rowData(IntensityExperiment)$ProteinId
  long <- melt(wide, id.vars = c("ProteinId"), variable.name = "IntensityColumn", value.name = "Intensity")
  long <- merge(long, colData(IntensityExperiment), by =  "IntensityColumn")
  as.data.table(long)
}


#' This function performs the log2 conversion and writes the imputed column
#' @param IntensityExperiment Output from constructSummarizedExperiment
#' @export prepare_prot_int
#' @import data.table

prepare_prot_int <- function(IntensityExperiment){
  prot.int <- get_long_protein_intensity(IntensityExperiment)
  stopifnot(dim(prot.int)[1]>0)
  prot.int <- as.data.table(prot.int)
  prot.int[, Imputed := 0L]
  prot.int[Intensity == 0, Imputed := 1L]
  prot.int[, log2NInt := 0.0]
  prot.int[Intensity > 0 , log2NInt := log2(Intensity)]
  
  prot.int 
}


