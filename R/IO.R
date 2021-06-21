#' This function writes a protein viz object which can be passed into MD's rails
#' database
#' 
#' @param outputFolder Path to output folder
#' @param IntensityExperiment Output from runGenericDiscovery
#' @export writeProteinViz
#' @importFrom data.table as.data.table
#' @importFrom SummarizedExperiment rowData
#' @importFrom jsonlite write_json unbox

writeProteinViz <- function(outputFolder, IntensityExperiment){
  
  comparisonStrings <- get_comparison_strings(rowData(IntensityExperiment))
  conditions <- unique(IntensityExperiment$Condition)
  
  proteinViz = list()
  
  for (comparison in comparisonStrings){
    
    # print(comparison)
    
    #get up and down
    condition1 <- get_condition_string(conditions, comparison, position = 1)
    condition2 <- get_condition_string(conditions, comparison, position = 2)
    
    # get statistics
    comparisonStatistics <- rowData(IntensityExperiment)
    comparisonStatistics <- filter_stats_table_on_comparison(comparisonStatistics, comparison)
    comparisonStatistics <- fill_out_missing_columns(comparisonStatistics)
    comparisonStatistics <- rename_comparison_statistics_export(comparisonStatistics)
    
    # get line
    sigProteinsPValue = comparisonStatistics[comparisonStatistics$AdjustedPValue<=0.05,]$PValue
    fdrLine <- max(sigProteinsPValue,na.rm = T)
    
    x = list(conditionComparison = unbox(comparison),
             up.condition = unbox(condition1),
             down.condition = unbox(condition2),
             fdrLimit =  unbox(fdrLine),
             data = as.data.table(comparisonStatistics))
    
    proteinViz[[comparison]] = x
  }
  
  names(proteinViz) <- NULL
  
  dir.create(outputFolder)
  path = file.path(outputFolder,"protein_viz.json")
  write_json(proteinViz,path, digits = NA, na = "null")
}





#' Use literal string match to find condition name in comparison string
#' @export detect_condition_string
detect_condition_string <- function(conditions, string){
  for (condition in conditions){
    if (grepl(condition, string, fixed = T)){
      return(condition)
    } 
  }
  stop("Couldn't match a condition the assay data and comparisons")
}

#' Use literal string match to find condition name in comparison string
#' @export get_condition_string
get_condition_string <- function(conditions, comparison, position = 1){
  
  stopifnot((position == 1) | (position == 2))
  
  match1 <- detect_condition_string(conditions,comparison) 
  comparison_remainder = gsub(match1, "",comparison, fixed = T)
  match2 <- detect_condition_string(conditions, comparison_remainder)   
  
  stopifnot(match2 != match1)
  
  position_match_1 <- gregexpr(pattern = match1,
                               comparison,
                               fixed = T)[[1]][1]
  
  position_match_2 <- gregexpr(pattern = match2,
                               comparison,
                               fixed = T)[[1]][1]
  if (position == 1){ # give the first match
    if (position_match_1 > position_match_2){
      return(match2)
    } else {
      return(match1)
    }
  } else if (position == 2) { # give the second match
    if (position_match_2 > position_match_1){
      return(match2)
    } else {
      return(match1)
    }
  }
  stop("Something went wrong, shouldn't be accessible. ")
}



#' This function takes the statistics table from limma and retrieves a single comparison's worth of statistics.
#' @export filter_stats_table_on_comparison
#' @import stringr
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
  
  statisticsTable <- statisticsTable[,c("ProteinId",statisticsColumns)]
  
  colnames(statisticsTable) <- gsub(str_c(" ",comparison), "",colnames(statisticsTable))
  
  return(statisticsTable)
}


#' This function ensures that all the missing columns are filled with NA's.
#' It should probably flag important columns that are missing in the future. 
#' @export fill_out_missing_columns
fill_out_missing_columns <- function(comparisonStatistics){
  
  columnsPossible=  c("ProteinId","GeneId","Description",
                       "logFC","adj.P.Val","P.Value",
                       "CI.L", "CI.R")
  
  columnsPresent = intersect(columnsPossible, colnames(comparisonStatistics))
  columnsAbsent <- columnsPossible[!(columnsPossible %in% columnsPresent)]
  comparisonStatistics[,columnsAbsent]<-NA
  
  comparisonStatistics
}


#' MD's protein viz needs specific attribute names so this function does that renaming
#' @export rename_comparison_statistics_export
#' @importFrom dplyr rename

rename_comparison_statistics_export <- function(comparisonStatistics){
  
  comparisonStatistics = as_tibble(comparisonStatistics) %>% 
    dplyr::rename(
      ProteinId = ProteinId,
      GeneName = GeneId,
      Description = Description,
      PValue = P.Value,
      AdjustedPValue = adj.P.Val,
      FoldChange = logFC,
      ConfLow = CI.L, 
      ConfHigh = CI.R
    )
  
  comparisonStatistics
} 


#' This functions works out which comparisons were calculated by limma.
#' Currently, it calculates all of them but it assumes order dependence so not trivial.
#' @export get_comparison_strings

get_comparison_strings <- function(results.quant){
  cols <- colnames(results.quant)
  cols <- cols[grepl("logFC ", colnames(results.quant))]
  comparison.strings <- gsub("logFC ", "", cols)
  comparison.strings
}
