
#' This function creates a list of intensity experiments pertaining to each
#' comparison calculated by the limma pipeline
#' @export listComparisonExperiments
#' @importFrom SummarizedExperiment rowData 
listComparisonExperiments <- function(CompleteIntensityExperiment){
  
  comparisons = 
    get_comparison_strings(rowData(CompleteIntensityExperiment))
  
  
  comparisonExperiments <- list()
  
  for (comparison in comparisons){
    comparisonExperiments[[comparison]] = createComparisonExperiment(
      CompleteIntensityExperiment, comparison
    )
  }
  
  comparisonExperiments
}

# filters out stats, assay and experiment design not column not 
# relevant to the specified comparison
#' @export filterComparisonColumns
filterComparisonColumns <- function(CompleteIntensityExperiment,
                                       comparison){
  
  conditions <- unique(CompleteIntensityExperiment$Condition)
  condition1 <- get_condition_string(conditions, comparison, 1)
  condition2 <- get_condition_string(conditions, comparison, 2)
  
  
  ### Deal with Statistics in rowData
  row.data <- rowData(CompleteIntensityExperiment)
  intensityColumns <- colnames(row.data)
  # we want the limma stats which have the comparison string
  # we want the replicates and imputed counts which have the conditions strings
  # we want the basic protein info cols (prot, gene, desc)
  cols.index = grepl("ProteinId", intensityColumns, fixed=T)+
    grepl("GeneId", intensityColumns, fixed=T)+
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
#' @export filterComparisonRows
filterComparisonRows <- function(comparisonIntensityExperiment){
  rowDataPossible <- rowData(comparisonIntensityExperiment)
  rowIndexValid <- !is.na(rowDataPossible$P.Value)
  rowIndexValid
  comparisonIntensityExperiment <-
    comparisonIntensityExperiment[rowIndexValid,]
  
  # How do we know if this breaks? 
  comparisonIntensityExperiment
}

#' This is a simple check to validate the intensity/cls vector match the statistics
#' @export validateComparison
validateComparison <- function(rse){
  assay.data <- dplyr::as_data_frame(assay(rse))
  row_data_fc <- rowData(rse)$FC
  calc_fc <- rowMeans(assay.data[,rse$GROUP == 0]) - rowMeans(assay.data[,rse$GROUP == 1])
  stopifnot((row_data_fc - calc_fc) < 2**(-10))
}

# This is an orchestrator for create each comparison table input for enrichment
# from a completed protein summarized experiment
#' @export createComparisonExperiment
createComparisonExperiment <- function(CompleteIntensityExperiment,
                                          comparison){
  
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
  conditions <- unique(comparisonIntensityExperiment$Condition)
  condition2 <- get_condition_string(conditions, comparison, 2)
  comparisonIntensityExperiment$GROUP = 
    comparisonIntensityExperiment$Condition == condition2
  
  
  validateComparison(comparisonIntensityExperiment)
  
  # need to set columns first
  rownames(comparisonIntensityExperiment) <-
    rowData(comparisonIntensityExperiment)$ProteinId
  
  ### Temporary Code while Enrichment Browser is in use ###
  #need to map column names correctly to what EnrichmentBrowser expects
  rowData(comparisonIntensityExperiment) = 
    dplyr::as_data_frame(rowData(comparisonIntensityExperiment)) %>%
    rename(FC = logFC,
           ADJ.PVAL = adj.P.Val)
  
  comparisonIntensityExperiment
}
