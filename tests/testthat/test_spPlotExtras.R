context("spPlotExtras")

##run test coloursFromTargets
test_that("check that coloursFromTargets outputs the expected result", {
  
  ###TEST1####
  #run function
  singlets <- getData(test_spCountsSng, "counts.cpm")
  output <- coloursFromTargets(col40(), singlets, rownames(singlets)[1:3])
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
})
