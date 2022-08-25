#' Create final results

#' @param IntensityExperiment output from `createSummarizedExperiment`
#' @param limmaStats Statistic table from limma output
#' @param longIntensityDT Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).
#' @param conditionComparisonMapping a dataframe matching different limma statistics comparison strings to up and down conditions

#' @export createResults

createResults <- function(IntensityExperiment, 
                         limmaStats,
                         longIntensityDT = longIntensityDT,
                         conditionComparisonMapping = conditionComparisonMapping){
  
  CompleteIntensityExperiment <- createCompleteIntensityExperiment(IntensityExperiment,
                                                                   limmaStats=limmaStats,
                                                                   longIntensityDT = longIntensityDT,
                                                                   conditionComparisonMapping = conditionComparisonMapping)
  
  list(IntensityExperiment=IntensityExperiment,
       CompleteIntensityExperiment=CompleteIntensityExperiment, 
       longIntensityDT=longIntensityDT)
}