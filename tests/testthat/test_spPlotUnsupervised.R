#context("spPlotUnsupervised")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[ ,s], testErcc[ ,s])
cObjMul <- spCounts(testCounts[ ,!s], testErcc[ ,!s])
uObj <- testUns
sObj <- testSwa

##run test plotUnsupervisedData
test_that("check that plotUnsupervisedData outputs the expected result", {
  
  ###TEST1####
  #run function
  output <- plotUnsupervisedData(uObj, cObjSng)
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
  expect_true(ncol(output) == 5)
  expect_true(nrow(output) == table(s)['TRUE'])
})

##run test plotUnsupervisedClass
test_that("check that plotUnsupervisedClass outputs the expected result", {
  
  #run function
  output <- plotUnsupervisedClass(uObj, cObjSng)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotUnsupervisedMarkers
test_that("check that plotUnsupervisedMarkers outputs the expected result", {
  
  #run function
  output <- plotUnsupervisedMarkers(uObj, cObjSng, markers = "a10")
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
