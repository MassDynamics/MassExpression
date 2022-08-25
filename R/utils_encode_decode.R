##########################
# Internal function to encode and decode condition names to avoi problems with makeContrats
##########################

#' Create sanitized condition names to be valid for input to makeContrasts
#' @param dt data.table with columns `Condition` to be encoded
#' @param condNames vector of at least one specifying all conditions names which should be encoded
#' @keywords internal
#' @noRd
condition_name_encoder <- function(dt, condNames) {
  dt_copy = copy(dt)
  conditionsDict <-list()
  
  for(cond in condNames){
    conditionDict <- data.table(original=dt_copy[, unique(get(cond))])
    set.seed(255)
    conditionDict[, safe := stringi::stri_rand_strings(.N, 5, pattern = "[A-Za-z]")]
    
    dt_copy[, original := get(cond)]
    dt_copy <- merge(dt_copy, conditionDict, by="original", all = T)
    dt_copy[, cond] = dt_copy$safe
    dt_copy[, `:=`(original = NULL, safe = NULL)]
    conditionsDict[[cond]] <- conditionDict
  }
  
  list(dt=dt_copy, conditionsDict=conditionsDict)
}


#' Decode condition names for intensity data to return initial names provided in input
#' @param dt data.table with columns `Condition` to be decoded
#' @param dict mapping dictionary to use for condition name decoding
#' @keywords internal
#' @noRd
condition_name_decode_intensity_data <- function(dt, dict){
  dt_copy = copy(dt)
  dt_copy <- merge(dt_copy, dict, by.x = "Condition", by.y = "safe", all.x = T)
  dt_copy[, `:=`(Condition = NULL)]
  setnames(dt_copy, "original", "Condition")
  dt_copy[, RunId := str_c(Condition, Replicate, sep = ".")]
  return(dt_copy)
}


#' Decode to initial condition names in limma table
#' @param dt limma table
#' @param dict mapping dictionary to use for condition name decoding
#' @keywords internal
#' @noRd
condition_name_decode_limma_table <- function(dt, dict){
  dt_copy <- copy(dt)
  dt_columns <- colnames(dt_copy)
  for (i in 1:dict[, .N]) {
    original <- dict[i, original]
    safe <- dict[i, safe]
    
    dt_columns <- str_replace(dt_columns, pattern = safe, replacement = original)
  }
  colnames(dt_copy) <- dt_columns
  return(dt_copy)
}

#' Decode to initial condition names in comparison mapping 
#' @param dt data.tabl with mapping of pairwise comparison performed
#' @param dict mapping dictionary to use for condition name decoding
#' @keywords internal
#' @noRd
condition_name_decode_comparison_mapping <- function(dt, dict, conditionSeparator){
  dt_copy <- copy(dt)
  for (i in 1:dict[, .N]) {
    original <- dict[i, original]
    safe <- dict[i, safe]
    
    dt_copy[, up.condition := str_replace(up.condition, pattern = safe, replacement = original)]
    dt_copy[, down.condition := str_replace(down.condition, pattern = safe, replacement = original)]
  }
  dt_copy[, comparison.string := paste(up.condition, down.condition, sep=conditionSeparator)]
  return(dt_copy)
}
