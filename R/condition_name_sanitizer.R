#' Create sanitized condition names to be valid for input to makeContrasts
#' @export

condition_name_encoder <- function(dt, condition_col_name) {
  dt_copy = copy(dt)
  dt_copy[, original := Condition]
  dt_copy[, Condition := make.names(original)]
  
  conditionsDict <- unique(dt_copy[,c("original", "Condition")])
  setnames(conditionsDict, "Condition", "safe")
  dt_copy[, `:=`(original = NULL)]
  
  list(dt=dt_copy, conditionsDict=conditionsDict)
}

#' Decode condition names to return initial names provided in input
#' @export
condition_name_decode_intensity_data <- function(dt, dict){
  dt_copy = copy(dt)
  dt_copy <- merge(dt_copy, dict, by.x = "Condition", by.y = "safe", all.x = T)
  dt_copy[, `:=`(Condition = NULL)]
  setnames(dt_copy, "original", "Condition")
  dt_copy[, RunId := str_c(Condition, Replicate, sep = ".")]
  return(dt_copy)
}


#' @export
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


#' @export
condition_name_decode_comparison_mapping <- function(dt, dict){
  dt_copy <- copy(dt)
  for (i in 1:dict[, .N]) {
    original <- dict[i, original]
    safe <- dict[i, safe]
    
    dt_copy[, up.condition := str_replace(up.condition, pattern = safe, replacement = original)]
    dt_copy[, down.condition := str_replace(down.condition, pattern = safe, replacement = original)]
  }
  dt_copy[, comparison.string := paste(up.condition, down.condition, sep="-")]
  return(dt_copy)
}
