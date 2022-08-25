library(MassExpression)
library(testthat)
library(data.table)

test_one_condition_encoding <- function(current, expected){
  test_that("special characters encoded correctly in long data.table",{
    result = min(expected$dt$Condition == current$dt$Condition)
    expect_true(as.logical(result))
  })

  test_that("mapping dictionary has one entry",{
    result = length(current$conditionsDict) == 1
    expect_true(as.logical(result))
  })
  
  test_that("original conditions are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Condition$original == current$conditionsDict$Condition$original)
    expect_true(as.logical(result))
  })
  
  test_that("safe conditions are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Condition$safe == current$conditionsDict$Condition$safe)
    expect_true(as.logical(result))
  })
}

test_one_condition_time_encoding <- function(current, expected){
  test_that("Time encoded correctly in long data.table",{
    result = min(expected$dt$Time == current$dt$Time)
    expect_true(as.logical(result))
  })
  
  test_that("mapping dictionary has one entry",{
    result = length(expected$conditionsDict) == 1
    expect_true(as.logical(result))
  })
  
  test_that("original conditions are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Time$original == current$conditionsDict$Time$original)
    expect_true(as.logical(result))
  })
  
  test_that("safe conditions are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Time$safe == current$conditionsDict$Time$safe)
    expect_true(as.logical(result))
  })
}

test_two_conditions_encoding <- function(current, expected){
  test_that("condition: 'Condition' is encoded correctly in long data.table",{
    result = min(expected$dt$Condition == current$dt$Condition)
    expect_true(as.logical(result))
  })
  
  test_that("condition: 'Time' is encoded correctly in long data.table",{
    result = min(expected$dt$Time == current$dt$Time)
    expect_true(as.logical(result))
  })
  
  test_that("mapping dictionary has two entries",{
    result = length(expected$conditionsDict) == 2
    expect_true(as.logical(result))
  })
  
  test_that("original conditions for 'Condition' are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Condition$original == current$conditionsDict$Condition$original)
    expect_true(as.logical(result))
  })
  
  test_that("safe conditions for 'Condition' are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Condition$safe == current$conditionsDict$Condition$safe)
    expect_true(as.logical(result))
  })
  
  test_that("original conditions for 'Time' are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Time$original == current$conditionsDict$Time$original)
    expect_true(as.logical(result))
  })
  
  test_that("safe conditions for 'Time' are correct in mapping dictionary",{
    result = min(expected$conditionsDict$Time$safe == current$conditionsDict$Time$safe)
    expect_true(as.logical(result))
  })
  
}



# Test
dt1 <- data.table(Condition = c(rep(">=1 meter",5), rep("<1 meter", 5)))
dt2 <- data.table(Condition = c(rep("≥1 meter",5), rep("<1 meter", 5)))
dt3 <- data.table(Condition = c(rep("$5># meter",5), rep("$10># meter", 5)))
dt4 <- data.table(Condition = c(rep("^>3 meter",5), rep("^>60 meter", 5)))
dt5 <- data.table(Condition = c(rep("3-5 meter",5), rep("5-10 meter", 5)))
dt6 <- data.table(Time = c(rep("2 hours",5), rep("5 hours", 5)))
dt_two_conditions <- data.table(Time = c(rep("2 hours",5), rep("5 hours", 5)),
                                Condition = c(rep("3-5 meter",5), rep("5-10 meter", 5)))

dt1_expected <- list(dt=data.table(Condition = c(rep("YJWrq",5), rep("bdcYa", 5))),
                     conditionsDict = list(Condition = 
                                             data.table(original = c(">=1 meter", "<1 meter"), 
                                                        safe = c("bdcYa", "YJWrq"))))

dt2_expected <- list(dt=data.table(Condition = c(rep("YJWrq",5), rep("bdcYa", 5))),
                     conditionsDict = list(Condition = 
                                             data.table(original = c("≥1 meter", "<1 meter"), 
                                                        safe = c("bdcYa", "YJWrq"))))

dt3_expected <- list(dt=data.table(Condition = c(rep("YJWrq",5), rep("bdcYa", 5))),
                     conditionsDict = list(Condition = 
                                             data.table(original = c("$5># meter", "$10># meter"), 
                                                        safe = c("bdcYa", "YJWrq"))))

dt4_expected <- list(dt=data.table(Condition = c(rep("bdcYa",5), rep("YJWrq", 5))),
                     conditionsDict = list(Condition = 
                                             data.table(original = c("^>3 meter", "^>60 meter"), 
                                                        safe = c("bdcYa", "YJWrq"))))

dt5_expected <- list(dt=data.table(Condition = c(rep("bdcYa",5), rep("YJWrq", 5))),
                     conditionsDict = list(Condition = 
                                             data.table(original = c("3-5 meter", "5-10 meter"), 
                                                        safe = c("bdcYa", "YJWrq"))))

dt6_expected <- list(dt=data.table(Time = c(rep("bdcYa",5), rep("YJWrq", 5))),
                     conditionsDict = list(Time = data.table(original = c("2 hours", "5 hours"), 
                                                             safe = c("bdcYa", "YJWrq"))))

dt_two_conditions_expected <- list(dt=data.table(Time = c(rep("bdcYa",5), rep("YJWrq", 5)),
                                                 Condition = c(rep("bdcYa",5), rep("YJWrq", 5))),
                     conditionsDict = list(Time = data.table(original = c("2 hours", "5 hours"), 
                                                 safe = c("bdcYa", "YJWrq")),
                                           Condition = data.table(original = c("3-5 meter", "5-10 meter"), 
                                                             safe = c("bdcYa", "YJWrq"))))

# Current
dt1_current <- MassExpression:::condition_name_encoder(dt1, "Condition")
dt2_current <- MassExpression:::condition_name_encoder(dt2, "Condition")
dt3_current <- MassExpression:::condition_name_encoder(dt3, "Condition")
dt4_current <- MassExpression:::condition_name_encoder(dt4, "Condition")
dt5_current <- MassExpression:::condition_name_encoder(dt5, "Condition")
dt6_current <- MassExpression:::condition_name_encoder(dt6, "Time")
dt_two_conditions_current <- MassExpression:::condition_name_encoder(dt_two_conditions, c("Time","Condition"))

## test condition length and encoding
test_one_condition_encoding(current = dt1_current, expected = dt1_expected)
test_one_condition_encoding(current = dt2_current, expected = dt2_expected)
test_one_condition_encoding(current = dt3_current, expected = dt3_expected)
test_one_condition_encoding(current = dt4_current, expected = dt4_expected)
test_one_condition_encoding(current = dt5_current, expected = dt5_expected)
test_one_condition_time_encoding(current = dt6_current, expected = dt6_expected)
test_two_conditions_encoding(current = dt_two_conditions_current, expected = dt_two_conditions_expected)

