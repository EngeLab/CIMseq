context("spPlotCounts")

##run test plotCountsData
test_that("check that plotCountsData outputs the expected result", {
  #This is largley tested in the estimateCells function
  
  ###TEST1####
  nc1 <- ncol(getData(test_spCountsSng, "counts"))
  nc2 <- ncol(getData(test_spCountsMul, "counts"))
  nc <- nc1 + nc2
  
  #run function
  output <- plotCountsData(test_spCountsSng, test_spCountsMul)
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
  expect_true(ncol(output) == 4)
  expect_true(nrow(output) == nc)
})

##run test plotCountsERCC
test_that("check that plotCountsERCC outputs the expected result", {
  
  #run function
  output <- plotCountsMarkers(test_spCountsSng, test_spCountsMul)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotCountsMarkers
test_that("check that plotCountsMarkers outputs the expected result", {
    
  #run function
  markers <- c("RNR2", "NME7")
  output <- plotCountsMarkers(test_spCountsSng, test_spCountsMul, markers)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
