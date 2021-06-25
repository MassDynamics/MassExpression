#' This function creates a summarized experiment object from a protein intensity table and experiment design. 
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @return A SummarizedExperiment object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList


constructSummarizedExperiment <- function(experimentDesign, proteinIntensities){
  
  stopifnot("Condition" %in% colnames(experimentDesign))
  stopifnot("Replicate" %in% colnames(experimentDesign))
  stopifnot("IntensityColumn" %in% colnames(experimentDesign))
  
  stopifnot(all(experimentDesign$IntensityColumn %in% colnames(proteinIntensities)))
  stopifnot("ProteinId" %in% colnames(proteinIntensities))
  
  
  rowDataPossible <-  c("ProteinId","GeneId","Description")
  rowDataPresent <- intersect(rowDataPossible, colnames(proteinIntensities))
  rowDataAbsent <- rowDataPossible[!(rowDataPossible %in% rowDataPresent)]
  
  assayData <- as.matrix(proteinIntensities[,experimentDesign$IntensityColumn])
  colnames(assayData) <- experimentDesign$IntensityColumn
  rownames(assayData) <- proteinIntensities$ProteinId
  
  rowFeatures <- proteinIntensities[,rowDataPresent,drop=FALSE] 
  rowFeatures[,rowDataAbsent] <- NA
  
  # construct summarized experiment object
  IntensityExperiment <- SummarizedExperiment(rowData = rowFeatures,
                                              assays= SimpleList(raw=assayData),
                                              colData = experimentDesign)
  
  colnames(IntensityExperiment) <- experimentDesign$IntensityColumn
  rownames(IntensityExperiment) <- rowData(IntensityExperiment)$ProteinId
  stopifnot(colnames(IntensityExperiment) == 
              SummarizedExperiment::colData(IntensityExperiment)$IntensityColumn)
  
  return(IntensityExperiment)
}