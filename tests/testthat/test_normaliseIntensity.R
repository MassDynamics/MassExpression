test_normalisation_median <- function(current, expected){
  test_that("not imputed median normalised intensities are as expected",{
    result = min(expected[!is.na(expected)] == current[!is.na(current)])
    expect_true(as.logical(result))
  })
  
  test_that("imputed median normalised intensities are as expected",{
    result = min(is.na(expected) == is.na(current))
    expect_true(as.logical(result))
  })
}

test_normalisation_none <- function(current, expected){
  test_that("not imputed not normalised intensities are as expected",{
    result = min(expected[!is.na(expected)] == current[!is.na(current)])
    expect_true(as.logical(result))
  })
  
  test_that("imputed not normalised intensities are as expected",{
    result = min(is.na(expected) == is.na(current))
    expect_true(as.logical(result))
  })
}

# Run tests
# Fake long intensities DT - 4 proteins, 2 replicates in each of 2 conditions
longIntensityDT <- data.table(ProteinId = c("P1","P2","P3","P4", 
                                             "P1","P2","P3","P4" ,
                                             "P1","P2","P3", "P4",
                                             "P1","P2","P3","P4"),
                               Condition = c("A","A","A","A",
                                             "A","A","A","A",
                                            "B","B","B","B",
                                            "B","B","B","B"),
                              Replicate = c(1,1,1,1,
                                            2,2,2,2,
                                            1,1,1,1,
                                            2,2,2,2),
                              log2NInt = c(0.0,1.0,4.3,3.9,
                                           4.6, 0, 3.3,3.4, 
                                           0.0, 2.0,2.3,2.5,
                                           2.5,0.0,2.5,1.0),
                              Imputed = c(1L, 0L, 0L, 0L, 
                                          0L, 1L, 0L, 0L, 
                                          1L, 0L, 0L, 0L, 
                                          0L, 1L, 0L, 0L))
longIntensityDT$SampleName <- paste(longIntensityDT$Condition, longIntensityDT$Replicate, sep = "-")

# we could have log2 which wrae = 1if the raw intensity was 1 - but not worrying about that
median_conditionA_R1 <- median(c(1,4.3,3.9)) 
median_conditionA_R2 <- median(c(4.6, 3.3, 3.4)) 
median_conditionB_R1 <- median(c(2,2.3,2.5)) 
median_conditionB_R2 <- median(c(2.5,2.5,1)) 

# expected median normalised intensities
longIntensityDT$log2IntNorm_median_expected <- c(NA,1.0-median_conditionA_R1, 4.3-median_conditionA_R1, 3.9-median_conditionA_R1,
                                          4.6-median_conditionA_R2, NA,3.3-median_conditionA_R2,3.4-median_conditionA_R2,
                                          NA, 2.0-median_conditionB_R1,2.3-median_conditionB_R1,2.5-median_conditionB_R1,
                                          2.5-median_conditionB_R2, NA, 2.5-median_conditionB_R2, 1.0-median_conditionB_R2)

# If no normalisation is applied 
longIntensityDT$log2Int_expected <- longIntensityDT$log2NInt

# tests - be aware thar data.table objects don't behave like other data frames - it's similar to python objects
longIntensityDT_median <- data.frame(normaliseIntensity(longIntensityDT,normalisationMethod='Median'))
longIntensityDT_none <- data.frame(normaliseIntensity(longIntensityDT,normalisationMethod='None'))

test_normalisation_median(current = longIntensityDT_median$log2NIntNorm, 
                          expected = longIntensityDT_median$log2IntNorm_median_expected)


test_normalisation_none(current = longIntensityDT_none$log2NInt, 
                          expected = longIntensityDT_none$log2Int_expected)




