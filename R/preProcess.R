#' @export
#' 
#' @importFrom stringr str_c


preProcess <- function(IntensityExperiment, 
                       normalisationMethod){
  
  ## TODO add checks for metadata 
  
  print("Starting pre-processing")
  metadataInfo <- metadata(IntensityExperiment)
  
  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  
  # Create valid condition names
  if(metadataInfo$experimentType$conditionOnly){
    condName <- metadataInfo$experimentType$condition1Name
    encodedCondition <- condition_name_encoder(dt = longIntensityDT, condName = condName)
  } else if(metadataInfo$experimentType$conditionTimeOrDose){
    
  }
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