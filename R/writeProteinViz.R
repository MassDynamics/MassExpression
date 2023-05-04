#' This function writes a protein viz object which can be passed into MD's rails
#' database
#' 
#' @param IntensityExperiment Output from runGenericDiscovery
#' @param outputFolder Path to output folder
#' 
#' @export writeProteinViz
#' @importFrom data.table as.data.table
#' @importFrom SummarizedExperiment rowData
#' @importFrom jsonlite write_json unbox

writeProteinViz <- function(IntensityExperiment, outputFolder){
  
  conditionComparisonMapping = metadata(IntensityExperiment)$conditionComparisonMapping
  comparisons <- conditionComparisonMapping$comparison.string
  
  stopifnot(length(comparisons)>0)
  
  proteinViz = list()
  
  compAvail <- colnames(rowData(IntensityExperiment))[grep("P.Value",colnames(rowData(IntensityExperiment)))]
  compAvail <- gsub("P.Value ","",compAvail)
  for (comparison in comparisons){
    
    print(comparison)
    if(!(comparison %in% compAvail)) next
    
    #get up and down
    up.condition <- getUpCondition(conditionComparisonMapping, comparison)
    down.condition <- getDownCondition(conditionComparisonMapping, comparison)
    
    # get statistics
    comparisonStatistics <- rowData(IntensityExperiment)
    comparisonStatistics <- filter_stats_table_on_comparison(comparisonStatistics, comparison)
    comparisonStatistics <- fill_out_missing_columns(comparisonStatistics)
    comparisonStatistics <- rename_comparison_statistics_export(comparisonStatistics)
    
    # get line
    sigProteinsPValue = comparisonStatistics[comparisonStatistics$AdjustedPValue<=0.05,]$PValue
    fdrLine <- max(sigProteinsPValue,na.rm = T)
    
    x = list(conditionComparison = unbox(comparison),
             up.condition = unbox(up.condition),
             down.condition = unbox(down.condition),
             fdrLimit =  unbox(fdrLine),
             data = as.data.frame(comparisonStatistics))
    
    proteinViz[[comparison]] = x
  }
  
  names(proteinViz) <- NULL
  
  print(outputFolder)
  dir.create(outputFolder, showWarnings = FALSE)
  outPath = file.path(outputFolder,"protein_viz.json")
  write_json(proteinViz,outPath, digits = NA, na = "null")
}





#' This function takes the statistics table from limma and retrieves a single comparison's worth of statistics.
#' @import stringr
#' @keywords internal
#' @noRd
filter_stats_table_on_comparison <- function(statisticsTable, comparison){
  example_col = str_c("P.Value ",comparison)
  stopifnot(example_col %in% colnames(statisticsTable))
  statisticsTable <- 
    
    # filter rows on valid statistics  
    statisticsTable[!is.na(statisticsTable[[example_col]]),]
  
  # filter columns for ProteinId and condition
  
  statisticsColumns <- str_c(c("logFC ","adj.P.Val ",
                                "P.Value ",
                                "CI.L ", "CI.R "), 
                              comparison)
  
  stopifnot(statisticsColumns %in% colnames(statisticsTable))
  
  statisticsTable <- statisticsTable[,c("ProteinId", "GeneName", "Description",
                                        statisticsColumns)]
  
  colnames(statisticsTable) <- gsub(str_c(" ",comparison), "",colnames(statisticsTable), fixed = TRUE)
  
  return(statisticsTable)
}


#' This function ensures that all the missing columns are filled with NA's.
#' It should probably flag important columns that are missing in the future. 
#' @keywords internal
#' @noRd
fill_out_missing_columns <- function(comparisonStatistics){
  
  columnsPossible=  c("ProteinId","GeneName","Description",
                       "logFC","adj.P.Val","P.Value",
                       "CI.L", "CI.R")
  
  columnsPresent = intersect(columnsPossible, colnames(comparisonStatistics))
  columnsAbsent <- columnsPossible[!(columnsPossible %in% columnsPresent)]
  comparisonStatistics[,columnsAbsent]<-NA
  
  comparisonStatistics
}


#' MD's protein viz needs specific attribute names so this function does that renaming
#' @importFrom dplyr rename
#' @keywords internal
#' @noRd

rename_comparison_statistics_export <- function(comparisonStatistics){
  
  initial_colnames <- colnames(comparisonStatistics)
  comparisonStatistics <- as_tibble(comparisonStatistics)
  colnames(comparisonStatistics) <- initial_colnames
  
  comparisonStatistics <- comparisonStatistics %>% 
    dplyr::rename(
      ProteinId = ProteinId,
      GeneName = GeneName,
      ProteinDescription = Description,
      PValue = P.Value,
      AdjustedPValue = adj.P.Val,
      FoldChange = logFC,
      ConfLow = CI.L, 
      ConfHigh = CI.R
    )
  
  comparisonStatistics
} 

