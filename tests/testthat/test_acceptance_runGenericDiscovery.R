library(MassExpression)
library(testthat)

test_raw_output <- function(current, expected, tolerance=10**-5){
  test_that("column names of assay in raw are the same",{
    result = min(colnames(expected) == colnames(current))
    expect_true(as.logical(result))
  })
  
  test_that("column names of rowData raw are the same",{
    result = min(colnames(rowData(expected)) == colnames(rowData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("column names of colData raw are the same",{
    result = min(colnames(colData(expected)) == colnames(colData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("sample names in colData raw are the same",{
    current_sample_names <- colData(current)$IntensityColumn
    expected_sample_names <- colData(expected)$IntensityColumn
    result = min(expected_sample_names == current_sample_names)
    expect_true(as.logical(result))
  })
  
  test_that("Condition columns in colData raw are the same",{
    current_condition <- colData(current)$Condition
    expected_condition <- colData(expected)$Condition
    result = min(current_condition == expected_condition)
    expect_true(as.logical(result))
  })
  
  test_that("protein ids in raw are the same",{
    current_proteinids <- rowData(current)$ProteinId
    expected_proteinids <- rowData(expected)$ProteinId
    result = min(expected_proteinids == current_proteinids)
    expect_true(as.logical(result))
  })
  
  test_that("raw intensities are equal", {
    current_int_vector <- c(assay(current))
    expected_int_vector <- c(assay(expected))
    
    approx_same = all.equal(current_int_vector, expected_int_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
}

test_complete_output <- function(current, expected, tolerance=10**-5){
  test_that("column names of assay in complete are the same",{
    result = min(colnames(expected) == colnames(current))
    expect_true(as.logical(result))
  })
  
  test_that("column names of rowData complete are the same",{
    result = min(colnames(rowData(expected)) == colnames(rowData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("column names of colData complete are the same",{
    result = min(colnames(colData(expected)) == colnames(colData(current)))
    expect_true(as.logical(result))
  })
  
  test_that("sample names in colData complete are the same",{
    current_sample_names <- colData(current)$IntensityColumn
    expected_sample_names <- colData(expected)$IntensityColumn
    result = min(expected_sample_names == current_sample_names)
    expect_true(as.logical(result))
  })
  
  test_that("Condition columns in colData complete are the same",{
    current_condition <- colData(current)$Condition
    expected_condition <- colData(expected)$Condition
    result = min(expected_condition == current_condition)
    expect_true(as.logical(result))
  })
  
  test_that("protein ids in complete are the same",{
    current_proteinids <- rowData(current)$ProteinId
    expected_proteinids <- rowData(expected)$ProteinId
    result = min(expected_proteinids == current_proteinids)
    expect_true(as.logical(result))
  })
  
  test_that("complete intensities are equal", {
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


# Run tests
design <- mq_lfq_data$design
intensities <- mq_lfq_data$intensities

listIntensityExperiments <- runGenericDiscovery(experimentDesign = design, 
                                                proteinIntensities = intensities, 
                                                NormalisationMethod = "None")

currentCompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
currentIntensityExperiment <- listIntensityExperiments$IntensityExperiment

load("../data/mq_lfq_output.RData")

test_raw_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)


test_complete_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)


test_limma_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)




