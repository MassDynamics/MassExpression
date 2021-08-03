test_imputeLFQ <- function(current, expected_imputed, expected_not_imputed){
  test_that("imputed values are as expected",{
    result = min(round(expected_imputed,4) == round(current[Imputed==1]$log2Norm,4))
    expect_true(as.logical(result))
  })
  
  test_that("not imputed values are as expected",{
    result = min(expected_not_imputed == current[Imputed==0]$log2Norm)
    expect_true(as.logical(result))
  })
}


# test

