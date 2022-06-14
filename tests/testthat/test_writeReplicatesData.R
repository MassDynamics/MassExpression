library(jsonlite)


test_writeReplicatesData <- function(current, expected){
test_that("ProteinID is the same",{
  result = min(expected$ProteinId == current$ProteinId)
  expect_true(as.logical(result))
})

test_that("Colnames are the same",{
  result = min(colnames(expected) == colnames(current))
  expect_true(as.logical(result))
})

test_that("condition name is the same",{
  result = min(expected$conditions[[1]]$name == current$conditions[[1]]$name)
  expect_true(as.logical(result))
  
  result = min(expected$conditions[[2]]$name == current$conditions[[2]]$name)
  expect_true(as.logical(result))
})

test_that("numberOfReplicateCount is the same",{
  result = min(expected$conditions[[1]]$numberOfReplicateCount == current$conditions[[1]]$numberOfReplicateCount)
  expect_true(as.logical(result))
  
  result = min(expected$conditions[[2]]$numberOfReplicateCount == current$conditions[[2]]$numberOfReplicateCount)
  expect_true(as.logical(result))
})

test_that("precentageOfReplicates is the same",{
  result = min(expected$conditions[[1]]$precentageOfReplicates == current$conditions[[1]]$precentageOfReplicates)
  expect_true(as.logical(result))
  
  result = min(expected$conditions[[2]]$precentageOfReplicates == current$conditions[[2]]$precentageOfReplicates)
  expect_true(as.logical(result))
})

}


# test
SampleDT <- data.table(ProteinId = c(rep("Prot1", 10), rep("Prot2", 10)),
                       GeneName = NA, Description = NA,
                       log2NIntNorm = rep(1.5, 20),
                       Condition = c(rep(c("A", "B"), each=5), rep(c("A", "B"), each=5)), 
                       Replicate = c(1:5, 1:5, 1:5, 1:5), 
                       Imputed =  c(0,0,0,0,0,
                                    0,0,0,1,1,
                                    1,0,0,0,0,
                                    1,1,0,0,0),
                       SampleName = "sample name here")

writeReplicateData(SampleDT, "../data/current")
current <- read_json("../data/current/protein_counts_and_intensity.json", simplifyVector=TRUE)

expected <- read_json("../data/protein_counts_and_intensity.json", simplifyVector=TRUE)

test_writeReplicatesData(current = current, expected = expected)

