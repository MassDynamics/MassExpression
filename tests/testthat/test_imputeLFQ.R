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


# Test Data
f_imputePosition <- 1.8
f_imputeStDev <- 0.3
myQuantDT <- data.table(Imputed = c(0L, 1L, 1L, 0L, 0L),
                        log2Norm = c(23, 0, 0, 24, 100), 
                        SampleName = c("S1", "S1", "S2", "S2", "S3"),
                        ProteinId = c("P1","P2","P1","P2","P1"))

pid <- myQuantDT[order(ProteinId)]

pid_imp <- imputeLFQ(myQuantDT = pid, id_type = "ProteinId",
                     int_type = "log2Norm", 
                     f_imputePosition = f_imputePosition, 
                     f_imputeStDev=f_imputeStDev)

# Expected
final_impute_sd <- 13.25104
final_impute_pos <- -30.50623
set.seed(255)
expected_imputed <- rnorm(2, final_impute_pos, final_impute_sd)
expected_not_imputed <- c(23, 100, 24)

# Run tests
test_imputeLFQ(current = pid_imp,
               expected_imputed=expected_imputed, 
               expected_not_imputed=expected_not_imputed)
