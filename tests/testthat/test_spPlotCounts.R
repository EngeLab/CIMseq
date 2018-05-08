context("spPlotCounts")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[ ,s], testErcc[ ,s])
cObjMul <- spCounts(testCounts[ ,!s], testErcc[ ,!s])
uObj <- testUns
sObj <- testSwa

##run test plotCountsData
test_that("check that plotCountsData outputs the expected result", {
  #This is largley tested in the estimateCells function
  
  ###TEST1####
  #run function
  output <- plotCountsData(cObjSng, cObjMul)
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
  expect_true(ncol(output) == 4)
  expect_true(nrow(output) == ncol(testCounts))
})

##run test plotCountsERCC
test_that("check that plotCountsERCC outputs the expected result", {
  
  #run function
  output <- plotCountsMarkers(cObjSng, cObjMul)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotCountsMarkers
test_that("check that plotCountsMarkers outputs the expected result", {
    
  #run function
  markers <- c("a1", "a10")
  output <- plotCountsMarkers(cObjSng, cObjMul, markers)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
