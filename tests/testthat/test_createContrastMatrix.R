library(MassExpression)
library(testthat)
library(data.table)

test_create_all_pairwise_comparisons <- function(current, expected){
  test_that("all pairwise comparison are create correctly with 2 levels provided", {
    expect_true(all(current == expected))
  })
  
  test_that("all pairwise comparison are create correctly with >2 levels provided", {
    expect_true(all(current == expected))
  })
}


# tests

twoLevels <- c("A", "B")
fourLevels <- c("A", "B", "C", "D")

expectedTwolevels <- data.table(left = "A", right = "B")
expectedFourlevels <- data.table(left = c("A","A","A","B","B","C"), 
                                right = c("B","C","D", "C", "D", "D"))

currentTwoLevels <- createAllPairwiseComparisons(twoLevels)
currentFourLevels <- createAllPairwiseComparisons(fourLevels)

test_create_all_pairwise_comparisons(currentTwoLevels, expectedTwolevels)
test_create_all_pairwise_comparisons(currentFourLevels, expectedFourlevels)
