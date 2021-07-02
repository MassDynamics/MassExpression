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
  stopifnot("SampleName" %in% colnames(experimentDesign))
  
  stopifnot(all(experimentDesign$SampleName %in% colnames(proteinIntensities)))
  stopifnot("ProteinId" %in% colnames(proteinIntensities))
  
  # prepare coldata
  design = design %>% group_by(Condition) %>% mutate(Replicate = row_number())
  
  # prepare rowdata
  rowDataPossible <-  c("ProteinId","GeneId","Description")
  rowDataPresent <- intersect(rowDataPossible, colnames(proteinIntensities))
  rowDataAbsent <- rowDataPossible[!(rowDataPossible %in% rowDataPresent)]
  
  assayData <- as.matrix(proteinIntensities[,experimentDesign$SampleName])
  colnames(assayData) <- experimentDesign$SampleName

  rownames(assayData) <- proteinIntensities$ProteinId
  
  rowFeatures <- proteinIntensities[,rowDataPresent,drop=FALSE] 
  rowFeatures[,rowDataAbsent] <- NA
  
  # construct summarized experiment object
  IntensityExperiment <- SummarizedExperiment(rowData = rowFeatures,
                                              assays= SimpleList(raw=assayData),
                                              colData = experimentDesign)
  
  colnames(IntensityExperiment) <- experimentDesign$SampleName
  rownames(IntensityExperiment) <- rowData(IntensityExperiment)$ProteinId
  stopifnot(colnames(IntensityExperiment) == 
              SummarizedExperiment::colData(IntensityExperiment)$SampleName)
  
  return(IntensityExperiment)
}