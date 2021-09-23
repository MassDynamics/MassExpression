#' This function creates a summarized experiment object from a protein intensity table and experiment design. 
#' 
#' @param experimentDesign data.frame. Experiment design provided in input by the user. Required columsn are: `SampleName` and `Condition`.
#' @param proteinIntensities data.frame. Wide matrix of intensities. Rows are proteins and columns are SampleNames. Required column: `ProteinId`. 
#' @param listMetadata list of metadata: `Species`, `LabellingMethod`, `NormalisationAppliedToAssay`.
#' @return A SummarizedExperiment object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr group_by mutate row_number


createSummarizedExperiment <- function(experimentDesign, proteinIntensities, listMetadata){
  
  if(!("Condition" %in% colnames(experimentDesign))){
    stop("'Condition' column is not available in the experiment design.")
  }
  experimentDesign$Condition <- as.character(experimentDesign$Condition)
  
  if(!("SampleName" %in% colnames(experimentDesign))){
    stop("'SampleName' column is not available in the experiment design.")
  }
  experimentDesign$SampleName <- as.character(experimentDesign$SampleName)
  
  if(!(all(experimentDesign$SampleName %in% colnames(proteinIntensities)))){
    missing_samples <- experimentDesign$SampleName[!(experimentDesign$SampleName %in% colnames(proteinIntensities))] 
    stop(paste0("The following SampleNames in the experiment design are missing from the protein intensities table: ", 
         paste(missing_samples, collapse = ", ")))
  }
  
  if(!("ProteinId" %in% colnames(proteinIntensities))){
    stop("'ProteinId' column is not available in the protein intensities table.")
  }
  proteinIntensities$ProteinId <- as.character(proteinIntensities$ProteinId)
  
  len_levels_condition <- length(names(table(experimentDesign$Condition)))
  if(len_levels_condition < 2){
    stop("Condition in the experiment design should have at least 2 groups.")
  }
  
  if(!(listMetadata$NormalisationAppliedToAssay %in% c("None", "Median"))){
    stop("normalisationMethod should be one one of: 'Median' or `None`.")
  }
  
  # prepare colData with Replicate column 
  experimentDesign = tibble::as_tibble(experimentDesign) %>% 
    group_by(Condition) %>% 
    mutate(Replicate = row_number())
  
  # prepare rowdata
  rowDataPossible <-  c("ProteinId","GeneName","Description")
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