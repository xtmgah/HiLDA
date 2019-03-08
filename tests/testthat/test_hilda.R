test_that("correctly load the test data", {
  inputFile <- system.file("data/sampleG.rdata", package = "HiLDA")
  load(inputFile)

  expect_equal(as.character(class(G)), "MutationFeatureData")
  expect_equal(G@type, "independent")
  expect_equal(G@flankingBasesNum, 5)
})



test_that("running the global test and the local test", {
  inputFile <- system.file("data/sampleG.rdata", package = "HiLDA")
  load(inputFile)
  K <- 3
  Param <- pmsignature::getPMSignature(G, K = K)
  hilda_global <- hilda_bayesfactor(inputG = G, inputParam = Param, refGroup = seq(1,20,2), n.iter = 100)
  hilda_local <- hilda_test(inputG = G, inputParam = Param, refGroup = seq(1,20,2), n.iter = 100)

  expect_equal(class(hilda_global), "rjags")
  expect_equal(class(hilda_local), "rjags")
})
