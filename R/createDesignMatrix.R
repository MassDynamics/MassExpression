#' Create the model design matrix
#' 
#' @param experimentDesign data.frame. Experiment design 
#' @param splineModel logical. TRUE to prepare an experiment design with splines.
#' @param df_splines Number of degrees of freedom to use to fit the spline.

#' @export

createDesignMatrix <- function(experimentDesign, splineModel = FALSE, df_splines = 5){
  
  designTypeInfo <- detectDesignType(experimentDesign)
  
  # Design for condition only
  if(designTypeInfo$conditionOnly){
    condName <- designTypeInfo$condition1Name
    info(default_logger, "Condition only experiment. Condition: ", condName)
    
    if(splineModel){
      stopifnot(length(designTypeInfo$condition1Levels) >= 6)
      design <- makeDesignConditionOnlySpline(designTypeInfo$experimentDesign, 
                                              conditionName = condName, 
                                              df_splines = df_splines)
      model_type <- "spline"
      comparisonGroups <- NA
    } else {
      stopifnot(length(designTypeInfo$condition1Levels) >= 2)
      design <- makeDesignConditionOnly(designTypeInfo$experimentDesign, 
                                        conditionName = condName)
      model_type <- "category"
      comparisonGroups <- colnames(design)
    }
    
  } else if(hasConditionWithTimeDose(experimentDesign)) {
    condNames <- getConditionName(experimentDesign)
    design <- makeDesignInteraction(experimentDesign, 
                                    timeDoseName =  designTypeInfo$condition2Name,
                                    conditionName = designTypeInfo$condition1Name)
    model_type <- "category-interaction"
    comparisonGroups <- colnames(design)[grep(":", colnames(design))]
  } else {
    print("Not implemented yet model with 2 generic conditions")
    design <- NULL
  }
  
  designTypeInfo$modelDesign <- design
  designTypeInfo$modelType <- model_type
  designTypeInfo$conditionGroups <- comparisonGroups
  return(designTypeInfo)
  
}
