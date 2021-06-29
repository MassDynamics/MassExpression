# functions for dealing with condition names and comparison mapping coming out of limma
# hold over until we deal with this better in MassExpression

#' This is called in Limma to create a mapping for comparison strings to conditions
#' used later when writing output.
#' @export assembleComparisonConditionMapping
assembleComparisonConditionMapping <- function(conditionComparisonMapping, seperator = " - "){
  
  colnames(conditionComparisonMapping) = c("up.condition", "down.condition")
  conditionComparisonMapping$comparison.string = str_c(conditionComparisonMapping$up.condition,
                                                       seperator,
                                                       conditionComparisonMapping$down.condition)
  
  return(conditionComparisonMapping)
}

#' This function uses the condition comparison mapping to get the up condition
#' @export getUpCondition
getUpCondition <- function(conditionComparisonMapping, comparison.string){
  comparison.string <- as.character(comparison.string)
  relevant_comparison_index = conditionComparisonMapping$comparison.string == comparison.string
  stopifnot(sum(relevant_comparison_index)==1)
  conditionComparisonMapping$up.condition[relevant_comparison_index]
}

#' This function uses the condition comparison mapping to get the down condition
#' @export getDownCondition
getDownCondition <- function(conditionComparisonMapping, comparison.string){
  comparison.string <- as.character(comparison.string)
  relevant_comparison_index = conditionComparisonMapping$comparison.string == comparison.string
  stopifnot(sum(relevant_comparison_index)==1)
  conditionComparisonMapping$down.condition[relevant_comparison_index]
}