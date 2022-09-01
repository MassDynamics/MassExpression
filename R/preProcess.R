#' Prepare initial data for statistical modeling
#' @description In this steps condition are encoded to safe strings; intensities are normalised when required and missing
#' values imputed by MNAR
#' 
#' @param IntensityExperiment SummarizedExperiment object as returned by `importData`
#' 
#' @export
#' 
#' @importFrom stringr str_c


preProcess <- function(IntensityExperiment){
  
  ## TODO add checks for correct metadata structure and its presence  
  metadataInfo <- metadata(IntensityExperiment)
  longIntensityDT <- initialiseLongIntensityDT(IntensityExperiment)
  
  info(MassExpressionLogger(), "Encode conditions to create safe names for processing")
  conditionNames <- c(metadataInfo$experimentType$condition1Name, 
                      metadataInfo$experimentType$condition2Name)
  encodedCondition <- condition_name_encoder(dt = longIntensityDT, condNames = conditionNames)
  
  longIntensityDT <- encodedCondition$dt
  conditionsDict <- encodedCondition$conditionsDict
  
  # Create Median Normalized Measurements in each Condition/Replicate
  if(!is.null(metadataInfo$inputMetadata)){
    normalisationMethod <- metadataInfo$inputMetadata$NormalisationAppliedToAssay
    info(MassExpressionLogger(), "Normalise intensities")
    longIntensityDT <- normaliseIntensity(longIntensityDT=longIntensityDT,
                                          normalisationMethod=normalisationMethod)
    
    # Imputation
    info(MassExpressionLogger(), "Impute")
    longIntensityDT <- imputeLFQ(longIntensityDT, 
                                 id_type = "ProteinId", 
                                 int_type = "log2NIntNorm",
                                 f_imputeStDev = 0.3,
                                 f_imputePosition = 1.8)
  } else {
    info(MassExpressionLogger(), "Assuming log2Intensities already pre-processed.")
    setnames(longIntensityDT, "log2NInt","log2NIntNorm")
  }
  
  ##TODO return in SummarisedExperiment shape not just long format as well
  list(longIntensityDT = longIntensityDT, conditionsDict = conditionsDict)
  
}




