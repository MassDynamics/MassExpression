#' This function creates a summarized experiment object from a protein intensity table and experiment design. 
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @return A SummarizedExperiment object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr group_by mutate row_number


constructSummarizedExperiment <- function(experimentDesign, proteinIntensities, listMetadata){
  
  stopifnot("Condition" %in% colnames(experimentDesign))
  stopifnot("SampleName" %in% colnames(experimentDesign))
  
  stopifnot(all(experimentDesign$SampleName %in% colnames(proteinIntensities)))
  stopifnot("ProteinId" %in% colnames(proteinIntensities))
  
  # prepare colData with Replicate column 
  experimentDesign = tibble::as_tibble(experimentDesign) %>% 
    group_by(Condition) %>% 
    mutate(Replicate = row_number())
  
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
                                              colData = experimentDesign, 
                                              metadata = listMetadata)
  
  colnames(IntensityExperiment) <- experimentDesign$SampleName
  rownames(IntensityExperiment) <- rowData(IntensityExperiment)$ProteinId
  stopifnot(colnames(IntensityExperiment) == 
              SummarizedExperiment::colData(IntensityExperiment)$SampleName)
  
  return(IntensityExperiment)
}