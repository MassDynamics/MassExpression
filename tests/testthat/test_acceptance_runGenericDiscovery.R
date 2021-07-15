library(MassExpression)
library(testthat)

test_raw_output <- function(current, expected, tolerance=10**-5){
  test_that("rawSE: column names of assay in raw are the same",{
    result = min(colnames(expected) == colnames(current))
    expect_true(as.logical(result))
  })
  
  test_that("rawSE: column names of rowData raw are the same",{
    result = min(colnames(rowData(expected)) == colnames(rowData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("rawSE: column names of colData raw are the same",{
    result = min(colnames(colData(expected)) == colnames(colData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("rawSE: sample names in colData raw are the same",{
    current_sample_names <- colData(current)$IntensityColumn
    expected_sample_names <- colData(expected)$IntensityColumn
    result = min(expected_sample_names == current_sample_names)
    expect_true(as.logical(result))
  })
  
  test_that("rawSE: Condition columns in colData raw are the same",{
    current_condition <- colData(current)$Condition
    expected_condition <- colData(expected)$Condition
    result = min(current_condition == expected_condition)
    expect_true(as.logical(result))
  })
  
  test_that("rawSE: protein ids in raw are the same",{
    current_proteinids <- rowData(current)$ProteinId
    expected_proteinids <- rowData(expected)$ProteinId
    result = min(expected_proteinids == current_proteinids)
    expect_true(as.logical(result))
  })
  
  test_that("rawSE: raw intensities are equal", {
    current_int_vector <- c(assay(current))
    expected_int_vector <- c(assay(expected))
    
    approx_same = all.equal(current_int_vector, expected_int_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
}

test_complete_output <- function(current, expected, tolerance=10**-5){
  test_that("completeSE: column names of assay in complete are the same",{
    result = min(colnames(expected) == colnames(current))
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: column names of rowData complete are the same",{
    result = min(colnames(rowData(expected)) == colnames(rowData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: column names of colData complete are the same",{
    result = min(colnames(colData(expected)) == colnames(colData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: sample names in colData complete are the same",{
    current_sample_names <- colData(current)$IntensityColumn
    expected_sample_names <- colData(expected)$IntensityColumn
    result = min(expected_sample_names == current_sample_names)
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: Condition columns in colData complete are the same",{
    current_condition <- colData(current)$Condition
    expected_condition <- colData(expected)$Condition
    result = min(expected_condition == current_condition)
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: protein ids in complete are the same",{
    current_proteinids <- rowData(current)$ProteinId
    expected_proteinids <- rowData(expected)$ProteinId
    result = min(expected_proteinids == current_proteinids)
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: protein ids are the same as the row names of the summarised experiment object",{
    current_proteinids <- rowData(current)$ProteinId
    current_rownames <- row.names(current)
    result = min(current_rownames == current_proteinids)
    expect_true(as.logical(result))
  })
  
  test_that("completeSE: complete intensities are equal", {
    current_int_vector <- c(assay(current))
    expected_int_vector <- c(assay(expected))
    approx_same = all.equal(current_int_vector, expected_int_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
}

test_limma_output <- function(current, expected, tolerance=10**-5){
  test_that("logFCs are equal", {
    current_vector <- rowData(current)[,"logFC AZD8931_resistant_SKBR3_AZDRc-Parental_SKBR3"]
    expected_vector <- rowData(expected)[,"logFC AZD8931_resistant_SKBR3_AZDRc-Parental_SKBR3"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("PValues are equal", {
    current_vector <- rowData(current)[,"P.Value AZD8931_resistant_SKBR3_AZDRc-Parental_SKBR3"]
    expected_vector <- rowData(expected)[,"P.Value AZD8931_resistant_SKBR3_AZDRc-Parental_SKBR3"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("Adj PValues are equal", {
    current_vector <- rowData(current)[,"adj.P.Val AZD8931_resistant_SKBR3_AZDRc-Parental_SKBR3"]
    expected_vector <- rowData(expected)[,"adj.P.Val AZD8931_resistant_SKBR3_AZDRc-Parental_SKBR3"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("AveExpr are equal", {
    current_vector <- rowData(current)[,"AveExpr"]
    expected_vector <- rowData(expected)[,"AveExpr"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
}


test_comparisons_output <- function(complete_current, comparison_current){
  test_that("intensities in comparison summarised experiments (SE) are the same as the complete SE",{
    result = min(complete_current == comparison_current)
    expect_true(as.logical(result))
  })
}

test_concordance_maxquant_output <- function(current_diff_fc, expected_diff_fc, 
                                             current_diff_pval, expected_diff_pval,
                                             tolerance=10**-5){
  test_that("logFC are not more different than a tolerance value",{
    approx_same = all.equal(current_diff_fc, expected_diff_fc, tolerance = tolerance)
    expect_true(approx_same)
  })
  
  test_that("P.Values are not more different than a tolerance value",{
    approx_same = all.equal(current_diff_pval, expected_diff_pval, tolerance = tolerance)
    expect_true(approx_same)
  })
}


###########
# Run tests
###########
make_long_wide_df <- function(matrix_prot, new_int_col = "IntME"){
  int_prot <- data.frame(matrix_prot)
  int_prot$ProteinId <- rownames(int_prot)
  long_int_me <- int_prot %>% pivot_longer(cols = c(LFQ.intensity.1_hu_C1:LFQ.intensity.6_hu_P3), 
                                           values_to = new_int_col, names_to = "SampleName")
  long_int_me$SampleName <- gsub("LFQ.intensity.","", long_int_me$SampleName)
  long_int_me$SampleName <- tolower(long_int_me$SampleName)
  return(long_int_me)
}


design <- mq_lfq_data$design
intensities <- mq_lfq_data$intensities
species <- mq_lfq_data$parameters[mq_lfq_data$parameters$X1 == "Species",2]
normMethod <- mq_lfq_data$parameters[mq_lfq_data$parameters$X1 == "UseNormalisationMethod",2]
labMethod <- mq_lfq_data$parameters[mq_lfq_data$parameters$X1 == "LabellingMethod",2]

listIntensityExperiments <- runGenericDiscovery(experimentDesign = design, 
                                                proteinIntensities = intensities, 
                                                normalisationMethod = normMethod,
                                                species = species, 
                                                labellingMethod = labMethod)

# Output to be checked
currentCompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
currentIntensityExperiment <- listIntensityExperiments$IntensityExperiment
currentcomparisonExperiments <- 
  listComparisonExperiments(currentCompleteIntensityExperiment)[[1]]

currentCompleteIntensityExperiment_longdf <- make_long_wide_df(data.frame(assay(currentCompleteIntensityExperiment)),
                                                               new_int_col = "Int")
currentComparisonExperiments_longdf <- make_long_wide_df(data.frame(assay(currentcomparisonExperiments)),
                                                         new_int_col = "IntComp")

compare_me <- currentComparisonExperiments_longdf %>% left_join(currentCompleteIntensityExperiment_longdf)


load("../data/mq_lfq_output.Rdata")

test_raw_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)


test_complete_output(current = currentCompleteIntensityExperiment, 
                expected = expectedCompleteIntensityExperiment)


test_limma_output(current = currentCompleteIntensityExperiment, 
                expected = expectedCompleteIntensityExperiment)


test_comparisons_output(complete_current = compare_me$Int,
                        comparison_current = compare_me$IntComp)

# Compare with maxquant workflow
load("../data/HER2_maxquant_workflow.RData")
current_new_run <- data_bench_maxquant %>% left_join(as_tibble(rowData(currentcomparisonExperiments)))
current_diff_fc_maxquant <- current_new_run$FC - current_new_run$Disco_logFC.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3
current_diff_pval_maxquant <- current_new_run$P.Value - current_new_run$Disco_P.Value.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3

test_concordance_maxquant_output(current_diff_fc = current_diff_fc_maxquant, 
                                 expected_diff_fc = expected_diff_fc_maxquant, 
                                 current_diff_pval = current_diff_pval_maxquant, 
                                 expected_diff_pval = expected_diff_pval_maxquant)

