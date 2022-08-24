#' @export
#' 
#' @importFrom stringr str_c


preProcess <- function(IntensityExperiment, 
                       normalisationMethod){
  
  print("Starting pre-processing")
  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  
  # Create valid condition names
  encodedCondition <- condition_name_encoder(dt = longIntensityDT)
  longIntensityDT <- encodedCondition$dt
  conditionsDict <- encodedCondition$conditionsDict
  
  # Create Median Normalized Measurements in each Condition/Replicate
  longIntensityDT <- normaliseIntensity(longIntensityDT=longIntensityDT,
                                        normalisationMethod=normalisationMethod)
  
  # Imputation
  longIntensityDT <- imputeLFQ(longIntensityDT, 
                               id_type = "ProteinId", 
                               int_type = "log2NIntNorm",
                               f_imputeStDev = 0.3,
                               f_imputePosition = 1.8)
  
  
  # RunId will be unique to a row wheraes replicate may not
  longIntensityDT[, RunId := str_c(Condition, Replicate, sep = ".")]
}