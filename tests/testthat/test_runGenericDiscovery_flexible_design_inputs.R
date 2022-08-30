library(MassExpression)
library(testthat)


test_missing_list_for_custom_comparisons <- function(tryObject){
  test_that("Throws error if comparisonType = 'custom' and customComparisonsList is missing", {
    expect_true(class(tryObject) == "try-error")
  })

  test_that("Error message is correct if comparisonType = 'custom' and customComparisonsList is missing", {
    expect_true(attributes(tryObject)$condition$message == 
                  "No comparisons of choice provided. See ?runGenericDiscovery argument customComparisonsList.")
  })
}

test_wrong_input_level_for_custom_comparisons <- function(tryObject){
  test_that("Throws error if customComparisonsList has levels not in design", {
    expect_true(class(tryObject) == "try-error")
  })
  
  test_that("Error message is correct if customComparisonsList has levels not in design", {
    expect_true(attributes(tryObject)$condition$message == 
                  "Some rightLevels in customComparisonsList are missing from the design. 
           Missing levels: stage2")
  })
}

test_wrong_input_level_for_ordered_comparisons <- function(tryObject){
  test_that("Throws error if orderConditionsList has levels not in design", {
    expect_true(class(tryObject) == "try-error")
  })
  
  test_that("Error message is correct if orderConditionsList has levels not in design", {
    expect_true(attributes(tryObject)$condition$message == 
                  "Some levels in orderConditionsList/baselineConditionList are missing from the design. \n           Missing levels: ulcerAS")
  })
}


test_custom_comparison_stats_output <- function(current, expected){
  test_that("Colnames of rowData of final results only contain comparison of interest",{
    expect_equal(current, expected)
    
  })
}

# 1. comparisonType = "custom", no custom DF provided

########
design <- fragpipe_data$design
intensities <- fragpipe_data$intensities
species <- fragpipe_data$parameters[fragpipe_data$parameters$X1 == "Species",2]
normMethod <- fragpipe_data$parameters[fragpipe_data$parameters$X1 == "UseNormalisationMethod",2]
labMethod <- fragpipe_data$parameters[fragpipe_data$parameters$X1 == "LabellingMethod",2]

# Custom
tryCustom <- try(runGenericDiscovery(experimentDesign = design, 
                                                proteinIntensities = intensities, 
                                                normalisationMethod = normMethod,
                                                species = species, 
                                                labellingMethod = labMethod, 
                                                comparisonType = "custom"),
                                silent = TRUE)

test_missing_list_for_custom_comparisons(tryCustom)


# Custom wrong levels
customComparisons_test <- list(Condition = data.frame(left = "control", right = "stage2"))

tryCustom <- try(runGenericDiscovery(experimentDesign = design, 
                                     proteinIntensities = intensities, 
                                     normalisationMethod = normMethod,
                                     species = species, 
                                     labellingMethod = labMethod, 
                                     comparisonType = "custom",
                                     customComparisonsList = customComparisons_test),
                 silent = TRUE)

test_wrong_input_level_for_custom_comparisons(tryCustom)


# Custom comparison correct
customComparisons_test <- list(Condition = data.frame(left = "control", right = "recurrence"))

custom_correct <- runGenericDiscovery(experimentDesign = design, 
                                     proteinIntensities = intensities, 
                                     normalisationMethod = normMethod,
                                     species = species, 
                                     labellingMethod = labMethod, 
                                     comparisonType = "custom",
                                     customComparisonsList = customComparisons_test)
CompleteExperiment <- custom_correct$CompleteIntensityExperiment
cols_current <- colnames(rowData(CompleteExperiment))

cols_expected <- c("ProteinId", "Description", "GeneName", "AveExpr", "t", "adj.P.Val", 
                   "logFC control - recurrence", 
                   "CI.L control - recurrence", "CI.R control - recurrence", 
                   "P.Value control - recurrence", 
                   "adj.P.Val control - recurrence", 
                   "NImputed: control", "NImputed: recurrence", "NImputed: remission", "NImputed: ulcer", 
                   "NReplicates: control", "NReplicates: recurrence", "NReplicates: remission", "NReplicates: ulcer")
test_custom_comparison_stats_output(cols_current, cols_expected)



# Baseline comparison correct
custom_correct <- runGenericDiscovery(experimentDesign = design, 
                                      proteinIntensities = intensities, 
                                      normalisationMethod = normMethod,
                                      species = species, 
                                      labellingMethod = labMethod, 
                                      comparisonType = "oneVSall", 
                                      baselineConditionList = list(Condition = "control"))
CompleteExperiment <- custom_correct$CompleteIntensityExperiment
cols_current <- colnames(rowData(CompleteExperiment))

cols_expected <- c("ProteinId", "Description", "GeneName", "AveExpr","F", "adj.P.Val", 
                   "logFC recurrence - control", 
                   "CI.L recurrence - control", "CI.R recurrence - control", 
                   "P.Value recurrence - control", "adj.P.Val recurrence - control", 
                   
                   "logFC remission - control", 
                   "CI.L remission - control", "CI.R remission - control", 
                   "P.Value remission - control", "adj.P.Val remission - control", 
                   
                   "logFC ulcer - control", 
                   "CI.L ulcer - control", "CI.R ulcer - control", 
                   "P.Value ulcer - control", "adj.P.Val ulcer - control",
                   
                   "NImputed: control", "NImputed: recurrence", "NImputed: remission", "NImputed: ulcer", 
                   "NReplicates: control", "NReplicates: recurrence", "NReplicates: remission", "NReplicates: ulcer")
test_custom_comparison_stats_output(cols_current, cols_expected)

# Order comparison correct
custom_correct <- runGenericDiscovery(experimentDesign = design, 
                                      proteinIntensities = intensities, 
                                      normalisationMethod = normMethod,
                                      species = species, 
                                      labellingMethod = labMethod, 
                                      comparisonType = "oneVSall",
                                      orderConditionsList = list(Condition = c("control","recurrence", "remission", "ulcer")))
CompleteExperiment <- custom_correct$CompleteIntensityExperiment
cols_current <- colnames(rowData(CompleteExperiment))

cols_expected <- c("ProteinId", "Description", "GeneName", "AveExpr","F", "adj.P.Val", 
                   "logFC recurrence - control", 
                   "CI.L recurrence - control", "CI.R recurrence - control", 
                   "P.Value recurrence - control", "adj.P.Val recurrence - control", 
                   
                   "logFC remission - control", 
                   "CI.L remission - control", "CI.R remission - control", 
                   "P.Value remission - control", "adj.P.Val remission - control", 
                   
                   "logFC ulcer - control", 
                   "CI.L ulcer - control", "CI.R ulcer - control", 
                   "P.Value ulcer - control", "adj.P.Val ulcer - control",
                   
                   "NImputed: control", "NImputed: recurrence", "NImputed: remission", "NImputed: ulcer", 
                   "NReplicates: control", "NReplicates: recurrence", "NReplicates: remission", "NReplicates: ulcer")
test_custom_comparison_stats_output(cols_current, cols_expected)


# Order comparison missing levels
tryObject <- try(runGenericDiscovery(experimentDesign = design, 
                                      proteinIntensities = intensities, 
                                      normalisationMethod = normMethod,
                                      species = species, 
                                      labellingMethod = labMethod, 
                                      comparisonType = "oneVSall",
                                      orderConditionsList = list(Condition = c("control","recurrence", "remission", "ulcerAS"))),
                      silent = TRUE)

test_wrong_input_level_for_ordered_comparisons(tryObject)
