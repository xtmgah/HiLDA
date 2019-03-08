library(testthat)
library(HiLDA)

test_check("HiLDA")

test_that("correctly load the test data", {
  inputFile <- system.file("data/sampleG.rdata", package = "HiLDA")
  load(inputFile)

  expect_equal(as.character(class(G)), "MutationFeatureData")
  expect_equal(G@type, "independent")
  expect_equal(G@flankingBasesNum, 5)
})
