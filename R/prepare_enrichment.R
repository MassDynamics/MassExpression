#' This function creates a list of intensity experiments pertaining to each
#' comparison calculated by the limma pipeline
#' @param CompleteIntensityExperiment SummarizedExperiment containing output from the differential expession analysis. 
#' @export listComparisonExperiments
#' @importFrom SummarizedExperiment rowData 

listComparisonExperiments <- function(CompleteIntensityExperiment){
  
  conditionComparisonMapping <- metadata(CompleteIntensityExperiment)$conditionComparisonMapping
  comparisons <- conditionComparisonMapping$comparison.string
  
  comparisonExperiments <- list()
  
  for (comparison in comparisons){
    comparisonExperiments[[comparison]] <- createComparisonExperiment(
      CompleteIntensityExperiment, comparison
    )
  }
  
  comparisonExperiments
}

#' This is an orchestrator for create each comparison table input for enrichment
#' from a completed protein summarized experiment
#' @param CompleteIntensityExperiment SummarizedExperiment containing limma results
#' @param comparison str. Character defining contrast of interest to be extracted. 
#' @noRd
createComparisonExperiment <- function(CompleteIntensityExperiment,
                                       comparison){
  
  
  conditionComparisonMapping <- metadata(CompleteIntensityExperiment)$conditionComparisonMapping
  up.condition <- getUpCondition(conditionComparisonMapping, comparison)
  down.condition <- getDownCondition(conditionComparisonMapping, comparison)
  
  
  # 1. Filter for statistics in rowData that we care about
  comparisonIntensityExperiment <- filterComparisonColumns(
    CompleteIntensityExperiment,
    comparison
  )
  
  # 2. Filter for rows where we haven't imputed too much
  comparisonIntensityExperiment <- filterComparisonRows(
    comparisonIntensityExperiment
  )
  
  # 3. Add class/group vectors so we know how to do enrichment
  comparisonIntensityExperiment$GROUP <- comparisonIntensityExperiment$Condition == down.condition
  
  #4. need to set columns first
  # rownames(comparisonIntensityExperiment) <- rowData(comparisonIntensityExperiment)$ProteinId
  
  # 5. Validate FC before assigning
  validateComparison(comparisonIntensityExperiment, comparison)
  
  ### Temporary Code while Enrichment Browser is in use ###
  #need to map column names correctly to what EnrichmentBrowser expects
  rowData(comparisonIntensityExperiment) <- 
    dplyr::as_tibble(rowData(comparisonIntensityExperiment)) %>%
    rename(FC = logFC,
           ADJ.PVAL = adj.P.Val)
  
  #6. Tag metadata
  metadata(comparisonIntensityExperiment)$up.condition <- up.condition
  metadata(comparisonIntensityExperiment)$down.condition <- down.condition
  
  comparisonIntensityExperiment
}



#' filters out stats, assay and experiment design not column not 
#' relevant to the specified comparison
#' @noRd
filterComparisonColumns <- function(CompleteIntensityExperiment,
                                    comparison){
  
  conditionComparisonMapping <- metadata(CompleteIntensityExperiment)$conditionComparisonMapping
  condition1 <- getUpCondition(conditionComparisonMapping, comparison)
  condition2 <- getDownCondition(conditionComparisonMapping, comparison)
  
  ### Deal with Statistics in rowData
  row.data <- rowData(CompleteIntensityExperiment)
  intensityColumns <- colnames(row.data)
  
  # we want the limma stats which have the comparison string
  # we want the replicates and imputed counts which have the conditions strings
  # we want the basic protein info cols (prot, gene, desc)
  cols.index <- grepl("ProteinId", intensityColumns, fixed=T)+
    grepl("GeneName", intensityColumns, fixed=T)+
    grepl("Description", intensityColumns, fixed=T)+
    grepl(str_c("Imputed: ",condition1), intensityColumns, fixed = T)+
    grepl(str_c("Replicates: ",condition1), intensityColumns, fixed = T)+
    grepl(str_c("Imputed: ",condition2), intensityColumns, fixed = T)+
    grepl(str_c("Replicates: ",condition2), intensityColumns, fixed = T)+
    grepl(comparison, intensityColumns, fixed = T)
  
  required.cols <- intensityColumns[which(as.logical(cols.index))]
  
  row.data <- row.data[,required.cols]
  
  # gotta rename the stats columns now
  colnames(row.data) <-
    gsub(str_c(" ",comparison), "",
         colnames(row.data),
         fixed=T)
  
  rowData(CompleteIntensityExperiment) <- row.data 
  
  ### Filter colData Experiment Design
  index <- CompleteIntensityExperiment$Condition %in% c(condition1, condition2)
  CompleteIntensityExperiment <- CompleteIntensityExperiment[,index]
  
  CompleteIntensityExperiment
}

#' This function filter out rows in the assay that have na statistics 
#' usually caused by the missing value filter for protein 
#' @noRd
filterComparisonRows <- function(comparisonIntensityExperiment){
  rowDataPossible <- rowData(comparisonIntensityExperiment)
  rowIndexValid <- !is.na(rowDataPossible$P.Value)
  # rowIndexValid
  comparisonIntensityExperiment <-
    comparisonIntensityExperiment[rowIndexValid,]
  
  # How do we know if this breaks? 
  comparisonIntensityExperiment
}

#' This is a simple check to validate the intensity/cls vector match the statistics
#' @noRd
validateComparison <- function(rse, comparison){
  assay.data <- dplyr::as_tibble(assay(rse))
  row_data_fc <- rowData(rse)$logFC
  calc_fc <- rowMeans(assay.data[,rse$GROUP == 0]) - rowMeans(assay.data[,rse$GROUP == 1])
  max_difference <- as.logical(max(((row_data_fc - calc_fc) < 2**(-10))))
  if(is.null(row_data_fc) | !max_difference){
    stop(paste0("logFC for comparison ", comparison," not validated."))
  }
}

# functions for dealing with condition names and comparison mapping coming out of limma

#' This function uses the condition comparison mapping to get the up condition
#' @param conditionComparisonMapping List of SummarizedExperiments, one per pairwise condition
#' @param comparison.string string defining the contrast of interest to extract
#' @noRd
getUpCondition <- function(conditionComparisonMapping, comparison.string){
  comparison.string <- as.character(comparison.string)
  relevant_comparison_index = conditionComparisonMapping$comparison.string == comparison.string
  stopifnot(sum(relevant_comparison_index)==1)
  conditionComparisonMapping$up.condition[relevant_comparison_index]
}

#' This function uses the condition comparison mapping to get the down condition
#' @param conditionComparisonMapping List of SummarizedExperiments, one per pairwise condition
#' @param comparison.string string defining the contrast of interest to extract
#' @noRd
getDownCondition <- function(conditionComparisonMapping, comparison.string){
  comparison.string <- as.character(comparison.string)
  relevant_comparison_index = conditionComparisonMapping$comparison.string == comparison.string
  stopifnot(sum(relevant_comparison_index)==1)
  conditionComparisonMapping$down.condition[relevant_comparison_index]
}