library(MassExpression)
library(here)
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
    current_sample_names <- colData(current)$SampleName
    expected_sample_names <- colData(expected)$SampleName
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
    current_sample_names <- colData(current)$SampleName
    expected_sample_names <- colData(expected)$SampleName
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
    current_vector <- rowData(current)[,"logFC AZD8931_resistant_SKBR3_AZDRc - Parental_SKBR3"]
    expected_vector <- rowData(expected)[,"logFC AZD8931_resistant_SKBR3_AZDRc - Parental_SKBR3"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("PValues are equal", {
    current_vector <- rowData(current)[,"P.Value AZD8931_resistant_SKBR3_AZDRc - Parental_SKBR3"]
    expected_vector <- rowData(expected)[,"P.Value AZD8931_resistant_SKBR3_AZDRc - Parental_SKBR3"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("Adj PValues are equal", {
    current_vector <- rowData(current)[,"adj.P.Val AZD8931_resistant_SKBR3_AZDRc - Parental_SKBR3"]
    expected_vector <- rowData(expected)[,"adj.P.Val AZD8931_resistant_SKBR3_AZDRc - Parental_SKBR3"]
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

test_limma_output_integers <- function(current, expected, tolerance=10**-5){
  test_that("logFCs are equal", {
    current_vector <- rowData(current)[,"logFC 1 - 2"]
    expected_vector <- rowData(expected)[,"logFC 1 - 2"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("PValues are equal", {
    current_vector <- rowData(current)[,"P.Value 1 - 2"]
    expected_vector <- rowData(expected)[,"P.Value 1 - 2"]
    approx_same = all.equal(current_vector, expected_vector, tolerance = tolerance)
    expect_true(approx_same) # tolerate small differences
  })
  
  test_that("Adj PValues are equal", {
    current_vector <- rowData(current)[,"adj.P.Val 1 - 2"]
    expected_vector <- rowData(expected)[,"adj.P.Val 1 - 2"]
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

test_qc_reports_exist <- function(output_folder){
  test_that("Test that the html QC report exists", {
    expect_true(file.exists(file.path(output_folder, "QC_Report.html")))
  })
  
  test_that("Test that the PDF QC report exists", {
    expect_true(file.exists(file.path(output_folder, "pdf","QC_Report.pdf")))
  })
  
  test_that("That that all separate QC reports exist", {
    qc_names <- get_names_qc()
    reports_exists <- sapply(qc_names, function(qc) file.exists(file.path(output_folder, paste0("QC_", qc, ".html"))))
    if(min(reports_exists) != 1) print(reports_exists)
    expect_true(min(reports_exists) == 1)
  })
  
}

###########
# Run tests
###########
make_long_wide_df <- function(matrix_prot, new_int_col = "IntME", with_integers=FALSE){
  int_prot <- data.frame(matrix_prot)
  int_prot$ProteinId <- rownames(int_prot)
  
  if(!with_integers){
  long_int_me <- int_prot %>% pivot_longer(cols = c(LFQ.intensity.1_hu_C1:LFQ.intensity.6_hu_P3), 
                                           values_to = new_int_col, names_to = "SampleName")
  long_int_me$SampleName <- gsub("LFQ.intensity.","", long_int_me$SampleName)
  long_int_me$SampleName <- tolower(long_int_me$SampleName)
  }else{
    long_int_me <- int_prot %>% pivot_longer(cols = c("X1":"X6"), 
                                             values_to = new_int_col, names_to = "SampleName")
  }
  return(long_int_me)
}

print("###########################")
print("test with MaxQuant output")
print("###########################")
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



intensitiesNAs <- intensities
intensitiesNAs[1:3,"LFQ.intensity.1_hu_C1"] <- NA
listIntensityExperimentsWithNas <- runGenericDiscovery(experimentDesign = design, 
                                                proteinIntensities = intensitiesNAs, 
                                                normalisationMethod = normMethod,
                                                species = species, 
                                                labellingMethod = labMethod)


# # QC reports
# print("Generate QC report")
# output_folder <- file.path(here(), "data/HER2-test-output/")
# print(paste0("Running tests from here:", output_folder))
# dir.create(output_folder, recursive=TRUE)
# qc_report <- system.file("rmd","QC_report.Rmd", package = "MassExpression")
# print(qc_report)
# generate_qc_report(listIntensityExperiments, output_folder = output_folder)
# generate_qc_report(listIntensityExperiments, output_folder = output_folder, format = "pdf")
# # 
# print("Generate separate QC reports")
# generate_separate_qc_reports(listIntensityExperiments, output_folder = output_folder)


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

load("../data/mq_lfq_output.RData")

test_raw_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)


test_complete_output(current = currentCompleteIntensityExperiment, 
                expected = expectedCompleteIntensityExperiment)


test_limma_output(current = currentCompleteIntensityExperiment, 
                expected = expectedCompleteIntensityExperiment)


test_comparisons_output(complete_current = compare_me$Int,
                        comparison_current = compare_me$IntComp)

##
print("###########################")
print("Test with NAs in intensity matrix")
print("###########################")
currentCompleteIntensityExperiment <- listIntensityExperimentsWithNas$CompleteIntensityExperiment
currentIntensityExperiment <- listIntensityExperimentsWithNas$IntensityExperiment
currentcomparisonExperiments <- 
  listComparisonExperiments(currentCompleteIntensityExperiment)[[1]]

currentCompleteIntensityExperiment_longdf <- make_long_wide_df(data.frame(assay(currentCompleteIntensityExperiment)),
                                                               new_int_col = "Int")
currentComparisonExperiments_longdf <- make_long_wide_df(data.frame(assay(currentcomparisonExperiments)),
                                                         new_int_col = "IntComp")

compare_me <- currentComparisonExperiments_longdf %>% left_join(currentCompleteIntensityExperiment_longdf)

test_raw_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)


test_complete_output(current = currentCompleteIntensityExperiment, 
                     expected = expectedCompleteIntensityExperiment)


test_limma_output(current = currentCompleteIntensityExperiment, 
                  expected = expectedCompleteIntensityExperiment)


test_comparisons_output(complete_current = compare_me$Int,
                        comparison_current = compare_me$IntComp)


################################
# Compare with output directly from maxquant workflow
load("../data/HER2_maxquant_workflow.RData")
current_new_run <- data_bench_maxquant %>% left_join(as_tibble(rowData(currentcomparisonExperiments)))
current_diff_fc_maxquant <- current_new_run$FC - current_new_run$Disco_logFC.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3
current_diff_pval_maxquant <- current_new_run$P.Value - current_new_run$Disco_P.Value.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3

test_concordance_maxquant_output(current_diff_fc = current_diff_fc_maxquant, 
                                 expected_diff_fc = expected_diff_fc_maxquant, 
                                 current_diff_pval = current_diff_pval_maxquant, 
                                 expected_diff_pval = expected_diff_pval_maxquant)

###### INTEGERS
###########################################################
##### Tests using integer input for Condition and SampleName
print("######################################################")
print("Tests using integer input for Condition and SampleName")
print("######################################################")

design <- mq_lfq_data$design
intensities <- mq_lfq_data$intensities

design <- design %>% mutate(Condition = ifelse(Condition %in% "AZD8931_resistant_SKBR3_AZDRc",1,2),
                            SampleName = 1:6)
colnames(intensities) <- c(1:6, "ProteinId")

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
                                                               new_int_col = "Int", with_integers = TRUE)
currentComparisonExperiments_longdf <- make_long_wide_df(data.frame(assay(currentcomparisonExperiments)),
                                                         new_int_col = "IntComp", with_integers = TRUE)

compare_me <- currentComparisonExperiments_longdf %>% left_join(currentCompleteIntensityExperiment_longdf)

load("../data/mq_lfq_output_integers.RData")

test_raw_output(current = currentIntensityExperiment, 
                expected = expectedIntensityExperiment)


test_complete_output(current = currentCompleteIntensityExperiment, 
                     expected = expectedCompleteIntensityExperiment)


test_limma_output_integers(current = currentCompleteIntensityExperiment, 
                  expected = expectedCompleteIntensityExperiment)


test_comparisons_output(complete_current = compare_me$Int,
                        comparison_current = compare_me$IntComp)

################################
# Compare with output directly from maxquant workflow
load("../data/HER2_maxquant_workflow_integers.RData")
current_new_run <- data_bench_maxquant %>% left_join(as_tibble(rowData(currentcomparisonExperiments)))
current_diff_fc_maxquant <- current_new_run$FC - current_new_run$Disco_logFC.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3
current_diff_pval_maxquant <- current_new_run$P.Value - current_new_run$Disco_P.Value.AZD8931_resistant_SKBR3_AZDRc...Parental_SKBR3

test_concordance_maxquant_output(current_diff_fc = current_diff_fc_maxquant, 
                                 expected_diff_fc = expected_diff_fc_maxquant, 
                                 current_diff_pval = current_diff_pval_maxquant, 
                                 expected_diff_pval = expected_diff_pval_maxquant)


########################
## QC report creation
#######################


