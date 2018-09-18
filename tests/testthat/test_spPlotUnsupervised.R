context("spPlotUnsupervised")

##run test plotUnsupervisedData
test_that("check that plotUnsupervisedData outputs the expected result", {
  
  ###TEST1####
  #run function
  output <- plotUnsupervisedData(test_spUnsupervised, test_spCountsSng)
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
  expect_true(ncol(output) == 4)
  expect_true(nrow(output) == ncol(getData(test_spCountsSng, "counts")))
})

##run test plotUnsupervisedClass
test_that("check that plotUnsupervisedClass outputs the expected result", {
  
  #run function
  output <- plotUnsupervisedClass(test_spUnsupervised, test_spCountsSng)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotUnsupervisedMarkers
test_that("check that plotUnsupervisedMarkers outputs the expected result", {
  
  #run function
  output <- plotUnsupervisedMarkers(
    test_spUnsupervised,
    test_spCountsSng,
    markers = "ACTB"
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
  
  #run function
  output <- plotUnsupervisedMarkers(
    test_spUnsupervised,
    test_spCountsSng,
    markers = c("ACTB", "GAPDH")
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
