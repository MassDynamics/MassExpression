#' This function creates a summarized experiment object from a protein intensity table and experiment design. 
#' 
#' @param experimentDesign XX.
#' @param proteinIntensities XX
#' @return A SummarizedExperiment object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr group_by mutate row_number


createSummarizedExperiment <- function(experimentDesign, proteinIntensities, listMetadata){
  
  if(!("Condition" %in% colnames(experimentDesign))){
    stop("'Condition' column is not available in the experiment design.")
  }
  
  if(!("SampleName" %in% colnames(experimentDesign))){
    stop("'SampleName' column is not available in the experiment design.")
  }
  
  if(!(all(experimentDesign$SampleName %in% colnames(proteinIntensities)))){
    missing_samples <- experimentDesign$SampleName[!(experimentDesign$SampleName %in% colnames(proteinIntensities))] 
    stop(paste0("The following SampleNames in the experiment design are missing from the protein intensities table: ", 
         paste(missing_samples, collapse = ", ")))
  }
  
  if(!("ProteinId" %in% colnames(proteinIntensities))){
    stop("'ProteinId' column is not available in the protein intensities table.")
  }
  
  len_levels_condition <- length(names(table(experimentDesign$Condition)))
  if(len_levels_condition < 2){
    stop("Condition in the experiment design should have at least 2 groups.")
  }
  
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