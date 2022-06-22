library(MassExpression)
library(testthat)


test_sanitize_invisible_characters <- function(current, expected){
  test_that("Test that invisible characters are removed and UNICODE characters converted to empty spaces.", {
    expect_equal(current, expected)
  })
}

test_sanitize_dataframe <- function(current, expected){
  test_that("Test that UNICODE/invisible characters in dataframe are handled correctly.", {
    expect_equal(current$SampleName, expected$SampleName)
    expect_equal(current$Condition, expected$Condition)
  })
}

##################################
# Examples sanitize single strings
##################################
badString1 <- "one \u200Btwo\u200B three"
badString2 <- "one \U00A0two\U00A0 three"
badString3 <- "\u00A0one two three\u00A0"
badString4 <- "Protein\u00A01"
string5 <- "Proteinu00A01"
string6 <- "Proteinuu200B"


expectedGoodString1 <- charToRaw("one two three")
expectedGoodString2 <- charToRaw("one  two  three")
expectedGoodString3 <- charToRaw("one two three")
expectedGoodString4 <- charToRaw("Protein 1")
expectedString5 <- charToRaw(string5)
expectedString6 <- charToRaw(string6)

sanitized_string1 <- charToRaw(sanitize_strings(badString1))
sanitized_string2 <- charToRaw(sanitize_strings(badString2))
sanitized_string3 <- charToRaw(sanitize_strings(badString3))
sanitized_string4 <- charToRaw(sanitize_strings(badString4))
sanitized_string5 <- charToRaw(sanitize_strings(string5))
sanitized_string6 <- charToRaw(sanitize_strings(string6))


# Test single strings
print("Test sanitize strings.")
test_sanitize_invisible_characters(current = sanitized_string1, expected = expectedGoodString1)
test_sanitize_invisible_characters(current = sanitized_string2, expected = expectedGoodString2)
test_sanitize_invisible_characters(current = sanitized_string3, expected = expectedGoodString3)
test_sanitize_invisible_characters(current = sanitized_string4, expected = expectedGoodString4)
test_sanitize_invisible_characters(current = sanitized_string5, expected = expectedString5)
test_sanitize_invisible_characters(current = sanitized_string6, expected = expectedString6)

##################################
# Examples sanitize dataframes
##################################
# Test that each row is sanitized
# Test that trailing and leading white spaces are remoced
# Test that combinations of invisible and unicode char are sanitized correctly

exampleBadDesign1 <- data.frame(SampleName = c("Sample1\u200B", "Sample2\u200B", 
                                              "Sample3\u00A0", 
                                              "Sample3 rep\u00A0\u200B"), 
                               Condition = c("\u00A0C1","\u200BC2", "C1", "C2")) 

exampleGoodDesign1 <- data.frame(SampleName = c("Sample1", "Sample2", 
                                               "Sample3", 
                                               "Sample3 rep"), 
                                Condition = c("C1","C2", "C1", "C2")) 


exampleDesign2 <- data.frame(SampleName = c(1,2,3,4), 
                                Condition = c(1,2,3,4)) 


# design 1
print("Test design strings to sanitize")
sanitized_dataframe <- sanitize_strings_in_dataframe(exampleBadDesign)
test_sanitize_dataframe(current = sanitized_dataframe, expected = exampleGoodDesign)

# design 2
print("Test design with no strings to sanitize")
sanitized_dataframe <- sanitize_strings_in_dataframe(exampleDesign2)
test_sanitize_dataframe(current = sanitized_dataframe, expected = exampleDesign2)

