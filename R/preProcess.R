#' Prepare initial data for statistical modeling
#' @description In this steps condition are encoded to safe strings; intensities are normalised when required and missing
#' values imputed by MNAR
#' 
#' @param IntensityExperiment SummarizedExperiment object as returned by `importData`
#' @param normalisationMethod str. One of 'None' or 'Median'
#' 
#' @export
#' 
#' @importFrom stringr str_c


preProcess <- function(IntensityExperiment, 
                       metadataInfo,
                       normalisationMethod){
  
  ## TODO add checks for metadata 
  
  print("Starting pre-processing...")
  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  
  print("Encode conditions to create safe names for processing")
  conditionNames <- c(metadataInfo$experimentType$condition1Name, 
                      metadataInfo$experimentType$condition2Name)
  encodedCondition <- condition_name_encoder(dt = longIntensityDT, condNames = conditionNames)
  
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
  
  ##TODO return in SummarisedExperiment shape not just long format as well
  list(longIntensityDT = longIntensityDT, conditionsDict = conditionsDict)
  
}