
#' Create pairwise comparisons
#' @param comparisonType Type of model / pairwise comparison. One of, "all", "oneVSall", "custom", "spline". 
#' Default to "all"  
#' @param contrastLevels vector with all condition levels to be considered for building contrasts
#' @param orderCondition vector of ordered condition levels
#' @param baselineInpuLevel vector of baseline level
#' @param chosenComparisons data.frame of chosen comparisons

#' @export createPairwiseComparisons

# customComparisons <- data.frame(left = c("T1", "T2"), right = c("T0", "T3"))

createPairwiseComparisons <- function(comparisonType = "all",
                                      condLevels, 
                                      conditionsDict, 
                                      orderCondition = NULL,
                                      baselineInpuLevel = NULL,
                                      customComparisonsCond = NULL){
  
  
  print("Encode vector of condition levels")
  condLevelsEncode <- encodeConditionComparisonsVec(allCondLevels = condLevels, 
                                                    condToEncode = condLevels, 
                                                    conditionsDict = conditionsDict)
  
  print("Encode vector of ordered condition levels provided in input")
  orderConditionEncode <- encodeConditionComparisonsVec(allCondLevels = condLevels, 
                                                        condToEncode = orderCondition, 
                                                        conditionsDict = conditionsDict)
  
  print("Encode baseline condition level provided in input")
  if(!is.null(baselineInpuLevel) & length(baselineInpuLevel) > 1){
    warning(paste0("baselineInpuLevel contains more than one level: ", 
                   paste(baselineInpuLevel, collapse = ","),
                   ". Only the first one will be considered for modeling."))
    baselineInpuLevel <- baselineInpuLevel[1]
  }
  baselineInpuLevelEncode <- encodeConditionComparisonsVec(allCondLevels = condLevels, 
                                                           condToEncode = baselineInpuLevel, 
                                                           conditionsDict = conditionsDict)
  
  customComparisonsCondEncode <- encodeCustomComparisonsDF(allCondLevels = condLevels, 
                                                           DFToEncode = customComparisonsCond, 
                                                           conditionsDict = conditionsDict)
  
  
  pairwiseComp <- switch(comparisonType, 
                         "all" = createAllPairwiseComparisons(condLevelsEncode),
                         "oneVSall" = createOneVsBaselineFromInput(condLevelsEncode, 
                                                                   orderConditionEncode, 
                                                                   baselineInpuLevelEncode),
                         "custom" = createOneVsOnePairwiseComparisons(condLevelsEncode, 
                                                                      customComparisonsCondEncode))
  
  return(pairwiseComp)

}

#' Create one vs baseline depending on input

#' @keywords internal
#' @noRd
#' 
createOneVsBaselineFromInput <- function(contrastLevels, 
                                         orderCondition,
                                         baselineInpuLevel){
  if(!is.null(orderCondition)){
    contrastLevels <- orderCondition
    baselineLevel <- orderCondition[1]
    print(paste0("Order of condition: ", paste(orderCondition, collapse = ",")))
    print(paste0("Baseline of condition: ", baselineLevel))
  } else if (!is.null(baselineInpuLevel)){
    baselineLevel <- baselineInpuLevel
    print("Order of condition not provided")
    print(paste0("Baseline of condition: ", baselineLevel))
  } else if (length(grep("[^0-9]", contrastLevels)) > 0){
    print(paste0("No baseline level was provided in input, defaulting to default order"))
    numbersInLevels <- gsub("[^0-9]", "", contrastLevels)
    numbersInLevels <- as.numeric(numbersInLevels)
    contrastLevels <- contrastLevels[order(numbersInLevels)]
    baselineLevel <- contrastLevels[1]
    print(paste0("Automatically detected order of condition: ", paste(contrastLevels, collapse = ",")))
    print(paste0("Automatically detected baseline of condition: ", baselineLevel))
  } else {
    stop("No baseline or order provided in input and condition levels don't have enough information for 
          automatic order detection.")
  }

  pairwiseComp <- createOneVsAllPairwiseComparisons(contrastLevels, baselineLevel)
  
  return(pairwiseComp)
  
} 

#' Create all pairwise comparisons

#' @keywords internal
#' @noRd
#' @import data.table

createAllPairwiseComparisons <- function(contrastLevels){
  cominationMat <- combn(x = contrastLevels, 2)
  pairwiseComp <- data.frame(t(cominationMat))
  colnames(pairwiseComp) <- c("left", "right")
  return(as.data.table(pairwiseComp))
}

#' Create one vs all pairwise comparisons

#' @keywords internal
#' @noRd
#' @import data.table

createOneVsAllPairwiseComparisons <- function(contrastLevels, baselineLevel){
  otherLevels <- contrastLevels[contrastLevels != baselineLevel]
  pairwiseComp <- data.frame(left = otherLevels, right = baselineLevel)
  return(as.data.table(pairwiseComp))
}


#' Create one vs one pairwise comparison

#' @keywords internal
#' @noRd
#' @import data.table

createOneVsOnePairwiseComparisons <- function(contrastLevels, chosenComparisons){
  if(is.null(chosenComparisons)){
    stop("No comparisons of choice provided. See ?runGenericDiscovery argument customComparisonsList.")
  }
  
  leftLevels <- chosenComparisons$left
  rightLevels <- chosenComparisons$right
  
  if(sum(leftLevels %in% contrastLevels) != length(leftLevels)){
    print(paste0("Some leftLevels of chosenComparison are not in the condition levels. Missing levels: ",
          paste(leftLevels[!(leftLevels %in% contrastLevels)], collapse = ",")))
  }
  
  if(sum(rightLevels %in% contrastLevels) != length(rightLevels)){
    print(paste0("Some rightLevels of chosenComparison are not in the condition levels. Missing levels: ",
                 paste(rightLevels[!(rightLevels %in% contrastLevels)], collapse = ",")))
  }
  
  return(as.data.table(chosenComparisons))
  
}



#' Create contrast matrix with all pairwise comparisons

#' @keywords internal
#' @noRd

createContrastsMatrix <- function(pairwiseComp, designMat, conditionSeparator = "-"){
  
  myContrasts <- apply(pairwiseComp, 1, function(x){
    paste0("Condition",x[1], conditionSeparator, "Condition",x[2])
    })
  
  contrastMatrix <- eval(as.call(c(
    as.symbol("makeContrasts"),
    as.list(myContrasts),
    levels = list(designMat)
  )))
  return(contrastMatrix)   
}



##### DELETE maybe

#' Create contrast matrix with all pairwise comparisons

#' @keywords internal
#' @noRd

createAllPairwiseContrastsMatrix <- function(pairwiseComp, conditionSeparator = "-"){
  
  myContrasts <- apply(pairwiseComp, 1, function(x) paste0(x[1], conditionSeparator, x[2]))
  
  contrast.matrix <- eval(as.call(c(
    as.symbol("makeContrasts"),
    as.list(myContrasts),
    levels = list(contrastLevels)
  )))
  return(contrast.matrix)   
}

#' Create contrast matrix with all vs one pairwise comparison

#' @keywords internal
#' @noRd

createOneVsAllContrastMatrix <- function(pairwiseComp, baselineLevel){
  
  if(!(baselineLevel %in% contrastLevels)){
    stop("baselineLevel is not part of the contrast levels")
  }
  
  pairwiseComp <- createOneVsAllPairwiseComparisons(contrastLevels, baselineLevel)
  myContrasts <- apply(pairwiseComp, 1, function(x) paste0(x[1], conditionSeparator, x[2]))
  
  contrast.matrix <- eval(as.call(c(
    as.symbol("makeContrasts"),
    as.list(myContrasts),
    levels = list(contrastLevels)
  )))
  return(contrast.matrix)   
}

#' Create contrast matrix with one vs one pairwise comparison

#' @keywords internal
#' @noRd

createOneVsOneContrastMatrix <- function(contrastLevels, levelUp, levelDown){
  
  if(!(levelUp %in% contrastLevels) | !(levelDown %in% contrastLevels)){
    stop("levelUp and levelDown are not part of the contrast levels.")
  }
  
  pairwiseComp <- createOneVsOnePairwiseComparisons(contrastLevels, levelUp, levelDown)
  myContrasts <- apply(pairwiseComp, 1, function(x) paste0(x[1], conditionSeparator, x[2]))

  contrast.matrix <- eval(as.call(c(
    as.symbol("makeContrasts"),
    as.list(myContrasts),
    levels = list(contrastLevels)
  )))
  return(contrast.matrix)   
}


