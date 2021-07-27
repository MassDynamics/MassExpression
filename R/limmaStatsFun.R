#' This function performs the differential expression analysis with limma including all pairwise comparisons 
#' using the condition provided
#' @param ID_type str. Column name of feature ID, e.g. `ProteinId`
#' @param int_type str. Column name containing intensities to use for the differential expression analysis, e.g. `log2NIntNorm`
#' @param condition_col_name str. Column name for the condition of interest, e.g. `Condition`. Used to build all pairwise comparisons.  
#' @param run_id_col_name str. Column name for the unique identifier of condition and replicate, e.g. `RunId`
#' @param rep_col_name str. Replicate column name. e.g. `Replicate`
#' @param funDT data.table. long table containing intensities for all proteins and samples. 
#' @param pairwise.comp chr. XX In the future this will allow the possibility of including only some pairwise comparisons of interest. 
#' @param all.comparisons logical. TRUE to run differential expression analysis for all pairwise comparisons.  
#' @param returnDecideTestColumn logical. If TRUE the row data of the `CompleteIntensityExperiment` will contain the output from 
#' `limma::decideTests`. 

#' @export limmaStatsFun
#' @import limma
#' @import foreach
#' @import Biobase 
#' @importFrom stringr str_order str_c str_sort
#' @importFrom data.table rbindlist dcast.data.table as.data.table

limmaStatsFun <- function(ID_type,
                            int_type,
                            condition_col_name,
                            run_id_col_name,
                            rep_col_name,
                            funDT,
                            pairwise.comp = NULL,
                            all.comparisons = TRUE,
                          returnDecideTestColumn) {
  
  
  # Create all possible pairwise comparisons using Condition
  if (all.comparisons) {
     comination_mat <- combn(x = funDT[, unique(get(condition_col_name))], 2)
    pairwise.comp <- foreach (i = 1:ncol(comination_mat), .packages="foreach") %do% {
      return(data.table(
        left = comination_mat[1,i],
        right = comination_mat[2,i]
      ))
    }
    pairwise.comp <- rbindlist(pairwise.comp)
  }
  
  funDT <- funDT[str_order(get(run_id_col_name), numeric = T)]
  funDT[, ID := str_c("ID.", get(ID_type))]
  
  # Filter protein
  filterDT <- prepareForFiltering(funDT=funDT, 
                             condition_col_name=condition_col_name, 
                             run_id_col_name=run_id_col_name)
  isPresent <- filterDT[repPC >= 0.5, unique(ID)]
  funDT <- funDT[ID %in% isPresent]
  
  # Create wide matrix of counts
  fun.formula <- as.formula(str_c("ID ~ ", run_id_col_name))
  int_matrix <-
    dcast.data.table(funDT, fun.formula, value.var = int_type)
  int_matrix <- int_matrix[str_order(ID, numeric = T)]
  intrawnames <- int_matrix[, ID]
  intcolnames <- colnames(int_matrix[, 2:ncol(int_matrix)])
  intcolnames <- str_sort(intcolnames, numeric = TRUE)
  
  int_matrix <- as.matrix(int_matrix[, intcolnames, with = FALSE])
  rownames(int_matrix) <- intrawnames
  
  ### Imputed value mask matrix ------
  imputed_matrix <-
    dcast.data.table(funDT, fun.formula, value.var = "Imputed")
  imputed_matrix <- imputed_matrix[str_order(ID, numeric = T)]
  imp_mat_rawnames <- imputed_matrix[, ID]
  imp_mat_colnames <- colnames(imputed_matrix[, 2:ncol(imputed_matrix)])
  imp_mat_colnames <- str_sort(imp_mat_colnames, numeric = TRUE)
  
  imputed_matrix <- as.matrix(imputed_matrix[, imp_mat_colnames, with = FALSE])
  rownames(imputed_matrix) <- imp_mat_rawnames
  isZ <- imputed_matrix == 1
  
  # Create sample information dataframe
  sample_dt <- as.data.frame(unique(funDT[, .(
    run_id = get(run_id_col_name),
    condition = get(condition_col_name),
    Replicate = get(rep_col_name)
  )]))
  #sample_dt$condition <- make.names(sample_dt$condition)
  rownames(sample_dt) <- sample_dt$run_id
  
  # Create features information dataframe
  feature_dt <- data.frame(ID = rownames(int_matrix), otherinfo = NA)
  
  rownames(feature_dt) <- feature_dt$ID
  
  # Create ExpressionSet object
  eset <- ExpressionSet(
    assayData = int_matrix,
    phenoData = AnnotatedDataFrame(sample_dt),
    featureData = AnnotatedDataFrame(feature_dt)
  )
  
  # DE model
  design.mat <- model.matrix(~ 0 + condition,
                             data = pData(eset))
  myContrasts = NULL
  
  condition_seperator = "-"
  
  for (irow in 1:pairwise.comp[, .N]) {
    left <- pairwise.comp[irow, left]
    right <- pairwise.comp[irow, right]
    
    newContrast <- str_c("condition",left, condition_seperator, "condition",right)
    myContrasts = c(myContrasts, newContrast)
  }
  
  
  conditionComparisonMapping <- assembleComparisonConditionMapping(
    pairwise.comp, 
    seperator = condition_seperator
    )
  
  
  contrast.matrix <- eval(as.call(c(
    as.symbol("makeContrasts"),
    as.list(myContrasts),
    levels = list(design.mat)
  )))
  
  # Fit linear models
  fit <- lmFit(eset, design.mat)
  fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
  
  stats <-
    topTable(fit2,
             number = nrow(fit2),
             sort.by = "none"
    )
  stats <- as.data.table(stats)
  stats <- stats[, .(ID, AveExpr, F, adj.P.Val)]
  
  one_model_stats <- extractOneModelStats(fitObject=fit2, 
                                           statsANOVA=stats,
                                           pairwiseComp=pairwise.comp,
                                           myContrasts=myContrasts,
                                          returnDecideTestColumn=returnDecideTestColumn)
  
  sep_models_stats <- fitSeparateModels(stats=stats, eset=eset, 
                                         pairwise.comp=pairwise.comp,
                                         funDT=funDT, 
                                         filterDT=filterDT,
                                         condition_col_name=condition_col_name, 
                                         run_id_col_name=run_id_col_name, 
                                        returnDecideTestColumn=returnDecideTestColumn)


  one_model_stats[, ID := str_replace_all(ID, "ID.", "")]
  setnames(one_model_stats, "ID", ID_type)
  
  sep_models_stats[, ID := str_replace_all(ID, "ID.", "")]
  setnames(sep_models_stats, "ID", ID_type)
  
  return(list(statsSepModels=sep_models_stats, 
              eset = eset, 
              conditionComparisonMapping = conditionComparisonMapping, 
              statsOneModel=one_model_stats))
}



#' @keywords internal
#' @noRd
extractOneModelStats <- function(fitObject, statsANOVA, pairwiseComp, myContrasts, returnDecideTestColumn){
  stats <- statsANOVA
  for (ipair in 1:pairwiseComp[, .N]) {
    left <- pairwiseComp[ipair, left]
    right <- pairwiseComp[ipair, right]
    myContrasts = str_c("condition",left, "-", "condition",right)
    s_dt <-
      as.data.table(topTable(
        fitObject,
        number = nrow(fitObject),
        sort.by = "p",
        confint = 0.95,
        coef = myContrasts
      ))
    s_dt <- s_dt[, .(ID, logFC, CI.L, CI.R, P.Value, adj.P.Val)]
    comparison <- str_replace_all(myContrasts, "condition", "")
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
fitSeparateModels <- function(statsANOVA, eset, pairwise.comp, funDT, 
                              filterDT, condition_col_name, run_id_col_name, 
                              returnDecideTestColumn){
  stats <- statsANOVA
  #### single pairwise comparisons
  for (ipair in 1:pairwise.comp[, .N]) {
    subsecting <- funDT[get(condition_col_name) %in% pairwise.comp[ipair, c(left, right)], unique(get(run_id_col_name))]
    
    # Filter for IDs that are not present in at least one experiment in pairwise manner
    isPresent <- filterDT[get(condition_col_name) %in% pairwise.comp[ipair, c(left, right)] & repPC >= 0.5, unique(ID)]
    
    if (length(isPresent) > 0) {
      eset_pair <- eset[rownames(eset) %in% isPresent, colnames(eset) %in% subsecting]
      
      
      
      design.mat <- model.matrix(~ 0 + condition,
                                 data = pData(eset_pair))
      
      
      myContrasts = NULL
      left <- pairwise.comp[ipair, left]
      right <- pairwise.comp[ipair, right]
      myContrasts = c(myContrasts,
                      str_c("condition",left, "-", "condition",right))
      contrast.matrix <- eval(as.call(c(
        as.symbol("makeContrasts"),
        as.list(myContrasts),
        levels = list(design.mat)
      )))
      
      fit <- lmFit(eset_pair, design.mat)
      fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
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
      
      comparison <- str_replace_all(myContrasts, "condition", "")
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