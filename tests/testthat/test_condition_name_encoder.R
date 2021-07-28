library(MassExpression)
library(testthat)
library(data.table)

test_name_encoding <- function(current, expected){
  test_that("special characters encoded correctly in long data.table",{
    result = min(expected$dt$Condition == current$dt$Condition)
    expect_true(as.logical(result))
  })

  test_that("original conditions are correct in mapping dictionary",{
    result = min(expected$conditionsDict$original == current$conditionsDict$original)
    expect_true(as.logical(result))
  })
  
  test_that("safe conditions are correct in mapping dictionary",{
    result = min(expected$conditionsDict$safe == current$conditionsDict$safe)
    expect_true(as.logical(result))
  })
}


# Test
dt1 <- data.table(Condition = c(rep(">=1 meter",5), rep("<1 meter", 5)))
dt2 <- data.table(Condition = c(rep("≥1 meter",5), rep("<1 meter", 5)))
dt3 <- data.table(Condition = c(rep("$5># meter",5), rep("$10># meter", 5)))
dt4 <- data.table(Condition = c(rep("^>3 meter",5), rep("^>60 meter", 5)))
dt5 <- data.table(Condition = c(rep("3-5 meter",5), rep("5-10 meter", 5)))

dt1_expected <- list(dt=data.table(Condition = c(rep("YJWrq",5), rep("bdcYa", 5))),
                     conditionsDict = data.table(original = c(">=1 meter", "<1 meter"), safe = c("bdcYa", "YJWrq")))

dt2_expected <- list(dt=data.table(Condition = c(rep("YJWrq",5), rep("bdcYa", 5))),
                     conditionsDict = data.table(original = c("≥1 meter", "<1 meter"), safe = c("bdcYa", "YJWrq")))

dt3_expected <- list(dt=data.table(Condition = c(rep("YJWrq",5), rep("bdcYa", 5))),
                     conditionsDict = data.table(original = c("$5># meter", "$10># meter"), safe = c("bdcYa", "YJWrq")))

dt4_expected <- list(dt=data.table(Condition = c(rep("bdcYa",5), rep("YJWrq", 5))),
                     conditionsDict = data.table(original = c("^>3 meter", "^>60 meter"), safe = c("bdcYa", "YJWrq")))

dt5_expected <- list(dt=data.table(Condition = c(rep("bdcYa",5), rep("YJWrq", 5))),
                     conditionsDict = data.table(original = c("3-5 meter", "5-10 meter"), safe = c("bdcYa", "YJWrq")))


# Current
dt1_current <- MassExpression:::condition_name_encoder(dt1)
dt2_current <- MassExpression:::condition_name_encoder(dt2)
dt3_current <- MassExpression:::condition_name_encoder(dt3)
dt4_current <- MassExpression:::condition_name_encoder(dt4)
dt5_current <- MassExpression:::condition_name_encoder(dt5)

## test condition length and encoding
test_name_encoding(current = dt1_current, expected = dt1_expected)
test_name_encoding(current = dt2_current, expected = dt2_expected)
test_name_encoding(current = dt3_current, expected = dt3_expected)
test_name_encoding(current = dt4_current, expected = dt4_expected)
test_name_encoding(current = dt5_current, expected = dt5_expected)

