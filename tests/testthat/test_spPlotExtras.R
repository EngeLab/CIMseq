#context("spPlotExtras")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[, s], testErcc[, s])
cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
uObj <- testUns
sObj <- testSwa

##run test coloursFromTargets
test_that("check that coloursFromTargets outputs the expected result", {
  
  ###TEST1####
  #run function
  output <- coloursFromTargets(col40(), testCounts[, s], c("a1", "a10", "a2"))
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
})
