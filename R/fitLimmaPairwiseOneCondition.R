#' This function performs the differential expression analysis with limma including all pairwise comparisons 
#' using the condition provided
#' @param featureIdType str. Column name of feature ID, e.g. `ProteinId`
#' @param intensityType str. Column name containing intensities to use for the differential expression analysis, e.g. `log2NIntNorm`
#' @param conditionColname str. Column name for the condition of interest, e.g. `Condition`. Used to build all pairwise comparisons.  
#' @param runIdColname str. Column name for the unique identifier of condition and replicate, e.g. `RunId`
#' @param repColname str. Replicate column name. e.g. `Replicate`
#' @param longIntensityDT data.table. long table containing intensities for all proteins and samples. 
#' @param pairwiseComparisons chr. XX In the future this will allow the possibility of including only some pairwise comparisons of interest. 
#' @param all.comparisons logical. TRUE to run differential expression analysis for all pairwise comparisons.  
#' @param returnDecideTestColumn logical. If TRUE the row data of the `CompleteIntensityExperiment` will contain the output from 
#' `limma::decideTests`. 
#' @param conditionSeparator string. String used to separate up and down condition in output. 

#' @export fitLimmaPairwiseOneCondition

#' @importFrom stringr str_order str_c str_sort
#' @importFrom data.table rbindlist dcast.data.table as.data.table

fitLimmaPairwiseOneCondition <- function(conditionColname,
                                     longIntensityDT,
                                     metadataExperiment,
                                     featureIdType,
                                      intensityType,
                                      runIdColname,
                                      repColname,
                                     fitSeparateModels,
                                      returnDecideTestColumn,
                                    conditionSeparator, 
                                      pairwiseComparisons) {
  
  # Reorder/create new protein id column for simplicity 
  # numeric = TRUE uses number as numbers and doesn't treat them as strings
  longIntensityDT <- longIntensityDT[str_order(get(runIdColname), numeric = T)]
  longIntensityDT[, ID := str_c(get(featureIdType))]
  
  print("Identify features with more than 50% missing values per condition.")
  # Tag for filtering 
  filterDT <- prepareForFiltering(longIntensityDT=longIntensityDT, 
                             conditionColname=conditionColname, 
                             runIdColname=runIdColname)
  isPresent <- filterDT[repPC >= 0.5, unique(ID)]
  longIntensityDT <- longIntensityDT[ID %in% isPresent]
  
  print("Create wide matrix of intensities/imputed")
  # From wide to long matrix of counts
  intensitiesMatrix <- MassExpression:::pivotDTLongToWide(longIntensityDT, 
                                         idCol = "ID",
                                         colNamesFrom = runIdColname, 
                                         fillValuesFrom = intensityType)
  
  imputedMatrix <- MassExpression:::pivotDTLongToWide(longIntensityDT, 
                                                          idCol = "ID",
                                                          colNamesFrom = runIdColname, 
                                                          fillValuesFrom = "Imputed")
  
  #if(!useImputed){
  #  intensitiesMatrix[imputedMatrix == 1] <- NA
  #}
  
  print("Create ExpressionSet object for limma")
  eset <- MassExpression:::createExpressionSetFromLongDT(longIntensityDT = longIntensityDT,
                                          intensitiesMatrix = intensitiesMatrix)


  print("Create design and contrasts matrices")
  designMat <- model.matrix(~ 0 + Condition,
                             data = pData(eset))
  
  contrastMatrix <- MassExpression:::createContrastsMatrix(pairwiseComp = pairwiseComparisons, 
                                       designMat = designMat, 
                                       conditionSeparator = conditionSeparator)
  
  conditionComparisonMapping <- assembleComparisonConditionMapping(
    pairwiseComparisons, 
    seperator = conditionSeparator
  )
  
  print("Fit linear models")
  modelFit <- fitLinearModelLimma(eset = eset, 
                                  designMat = designMat, 
                                  contrastMatrix = contrastMatrix, 
                                  metadataExperiment = metadataExperiment)
  
  fitObject <- modelFit$fitObject
  statsANOVA <- modelFit$statsANOVA

  if(fitSeparateModels){
    print("fitSeparateModels")
    resultsModel <- fitSeparateModels(statsANOVA=statsANOVA, eset=eset, 
                                      pairwiseComparisons=pairwiseComparisons,
                                      longIntensityDT=longIntensityDT, 
                                      filterDT=filterDT,
                                      conditionColname=conditionColname, 
                                      runIdColname=runIdColname, 
                                      returnDecideTestColumn=returnDecideTestColumn, 
                                      conditionSeparator=conditionSeparator)
  } else { 
    print("extractOneModelStats")
    resultsModel <- MassExpression:::extractOneModelStats(fitObject=fitObject, 
                                             statsANOVA=statsANOVA,
                                             pairwiseComparisons=pairwiseComparisons,
                                            returnDecideTestColumn=returnDecideTestColumn, 
                                            conditionSeparator=conditionSeparator)
  }

  setnames(resultsModel, "ID", featureIdType)
  
  return(list(resultsModel=resultsModel, 
              eset = eset,
              conditionComparisonMapping = conditionComparisonMapping))
}



#' @keywords internal
#' @noRd
extractOneModelStats <- function(fitObject, 
                                 statsANOVA, 
                                 pairwiseComparisons,
                                 returnDecideTestColumn, 
                                 conditionSeparator){
  stats <- statsANOVA
  for (ipair in 1:pairwiseComparisons[,.N]) {
    left <- pairwiseComparisons[ipair, left]
    right <- pairwiseComparisons[ipair, right]
    myContrasts = str_c("Condition",left, conditionSeparator, "Condition",right)
    s_dt <-
      as.data.table(topTable(
        fitObject,
        number = nrow(fitObject),
        sort.by = "p",
        confint = 0.95,
        coef = myContrasts
      ))
    s_dt <- s_dt[, .(ID, logFC, CI.L, CI.R, P.Value, adj.P.Val)]
    comparison <- str_replace_all(myContrasts, "Condition", "")
    colnames(s_dt)[colnames(s_dt) != "ID"] <-
      str_c(colnames(s_dt)[colnames(s_dt) != "ID"], comparison, sep = " ")
    stats <- merge(stats, s_dt, by = "ID", all = T)
  }
  
  if(returnDecideTestColumn){
    # Decide test
    dtest <- decideTests(fitObject)
    pid <- rownames(dtest)
    dtest <- as.data.table(dtest)
    colnames(dtest) <- paste0("decide_", colnames(dtest))
    dtest$ID <- pid
    stats <- merge(stats, dtest, by = "ID", all = T)
  }
  
  return(stats)
}


#' Fit separate models one for each parwise comparison
#' @keywords internal
#' @noRd
fitSeparateModels <- function(statsANOVA, eset, pairwiseComparisons, longIntensityDT, 
                              filterDT, conditionColname, runIdColname, 
                              returnDecideTestColumn, conditionSeparator){
  stats <- statsANOVA
  #### single pairwise comparisons
  for (ipair in 1:pairwiseComparisons[,.N]) {
    subsecting <- longIntensityDT[get(conditionColname) %in% pairwiseComparisons[ipair, c(left, right)], unique(get(runIdColname))]
    
    # Filter for IDs that are not present in at least one experiment in pairwise manner
    isPresent <- filterDT[get(conditionColname) %in% pairwiseComparisons[ipair, c(left, right)] & repPC >= 0.5, unique(ID)]
    
    if (length(isPresent) > 0) {
      eset_pair <- eset[rownames(eset) %in% isPresent, colnames(eset) %in% subsecting]
      
      design.mat <- model.matrix(~ 0 + Condition,
                                 data = pData(eset_pair))
      
      myContrasts = NULL
      left <- pairwiseComparisons[ipair, left]
      right <- pairwiseComparisons[ipair, right]
      myContrasts = c(myContrasts,
                      str_c("Condition",left, conditionSeparator, "Condition",right))
      contrastMatrix <- eval(as.call(c(
        as.symbol("makeContrasts"),
        as.list(myContrasts),
        levels = list(design.mat)
      )))
      
      fit <- lmFit(eset_pair, design.mat)
      fit2 <- contrasts.fit(fit, contrasts = contrastMatrix)
      fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
      
      s_dt <-
        as.data.table(topTable(
          fit2,
          number = nrow(fit2),
          sort.by = "p",
          confint = 0.95,
          coef = myContrasts
        ))
      s_dt <- s_dt[, .(ID, logFC, CI.L, CI.R, P.Value, adj.P.Val)]
      
      comparison <- str_replace_all(myContrasts, "Condition", "")
      colnames(s_dt)[colnames(s_dt) != "ID"] <-
        str_c(colnames(s_dt)[colnames(s_dt) != "ID"], comparison, sep = " ")
      stats <- merge(stats, s_dt, by = "ID", all = T)
      
      if(returnDecideTestColumn){
        dtest <- decideTests(fit2) 
        pid <- rownames(dtest)
        dtest <- as.data.table(dtest)
        colnames(dtest) <- paste0("decide_", colnames(dtest))
        dtest$ID <- pid
        stats <- merge(stats, dtest, by = "ID", all = T)
      }
      
    }
  }
  return(stats)
}


#' Create a mapping for comparison strings to conditions used later when writing output.
#' @param conditionComparisonMapping List of SummarizedExperiments, one per pairwise condition
#' @param separator character to use as up/down condition separator 
#' @noRd
#' @keywords internal
assembleComparisonConditionMapping <- function(conditionComparisonMapping, seperator = " - "){
  
  colnames(conditionComparisonMapping) = c("up.condition", "down.condition")
  conditionComparisonMapping$comparison.string = str_c(conditionComparisonMapping$up.condition,
                                                       seperator,
                                                       conditionComparisonMapping$down.condition)
  
  return(conditionComparisonMapping)
}



#' Fit linear models using limma

#' @export fitLinearModelLimma
#' @import limma
#' 
fitLinearModelLimma <- function(eset, designMat, contrastMatrix, metadataExperiment){
  
  repMeas <- metadataExperiment$experimentType$hasRepMeas
  techRepl <- metadataExperiment$experimentType$hasTechRepl
  
  consCorr <- NULL
  if(repMeas){
    repl <- pData(eset)$Subject
    corfit <- duplicateCorrelation(eset, block = repl)
    consCorr <- ifelse(corfit$consensus < 0, consCorr, corfit$consensus)
    print(paste0("Consensus correlation of repeated subject measures:", consCorr))
  } else if (techRepl){
    repl <- pData(eset)$TechRepl
    corfit <- duplicateCorrelation(eset, block = repl)
    consCorr <- ifelse(corfit$consensus < 0, consCorr, corfit$consensus)
    print(paste0("Consensus correlation of technical replicates:", consCorr))
  } else {
    print("No repeated measurements or technical replicates in experiment.")
  }
  
  if(!is.null(consCorr)){
    fit <- lmFit(eset, designMat, block = repl, cor = consCorr)
  }else{
    fit <- lmFit(eset, designMat)  
  }
  
  fit2 <- contrasts.fit(fit, contrasts = contrastMatrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
  
  statsANOVA <-
    topTable(fit2,
             number = nrow(fit2),
             sort.by = "none"
    )
  statsANOVA <- as.data.table(statsANOVA)
  
  stat_select <- ifelse(ncol(contrastMatrix) == 1, "t", "F")
  statsANOVA <- statsANOVA[, .(ID, AveExpr, get(stat_select), adj.P.Val)]
  colnames(statsANOVA) <- c("ID", "AveExpr", stat_select, "adj.P.Val")
  
  list(fitObject = fit2, statsANOVA = statsANOVA)
  
}

#' Extract colData given the longIntensityDT 

#' @noRd
#' @keywords internal
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' 

createExpressionSetFromLongDT <- function(longIntensityDT, intensitiesMatrix){
  allColumns <- colnames(longIntensityDT)
  colnamesInputDesign <- allColumns[allColumns %in% c("SampleName", "Condition", 
                                                      "Time", "Dose", "Subject", 
                                                      "TechRepl","RunId", "Replicate")]
  extractedColData <- as.data.frame(unique(longIntensityDT[,..colnamesInputDesign]))
  rownames(extractedColData) <- extractedColData$RunId
  
  # Create features information dataframe
  featureDF <- data.frame(ID = rownames(intensitiesMatrix), otherinfo = NA)
  rownames(featureDF) <- featureDF$ID
  
  # Create ExpressionSet object
  eset <- ExpressionSet(
    assayData = intensitiesMatrix,
    phenoData = AnnotatedDataFrame(extractedColData),
    featureData = AnnotatedDataFrame(featureDF)
  )
  
  return(eset)
}
