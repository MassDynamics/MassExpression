#' This function normalises (based on `normalisationMethod`) the log2 intensities stored in the `log2NInt` column initialised by `initialiseLongIntensityDT`. 
#' 
#' @param longIntensityDT Output from `initialiseLongIntensityDT`
#' @param normalisationMethod logical. Whether to perform median normalisation by condition. 
#' @export normaliseIntensity
#' @return log2NIntNorm If `normalisationMethod` is not `None` a new column is added containing the `log2NInt` normalised intensities 
#' @details Normalisation is only applied to non imputed intensities. `log2NIntNorm` is NA for imputed intensities. 

normaliseIntensity <- function(longIntensityDT, normalisationMethod){
  normaliseDT <- longIntensityDT
  # Create Median Normalized Measurements in each Condition/Replicate
  if (normalisationMethod == "Median"){
    normaliseDT[Imputed == F, log2NIntNorm := log2NInt - median(log2NInt), by = SampleName]
  }else if(normalisationMethod == "None") {
    normaliseDT[, log2NIntNorm := log2NInt]
  } else {
    stop(paste0("Normalisation method ",normalisationMethod," not availble.  Available methods are: `Median` and `None`"))
  }
  return(normaliseDT)
}
