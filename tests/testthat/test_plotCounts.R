context("plotCounts")

##run test plotCountsData
test_that("check that plotCountsData outputs the expected result", {
  #This is largley tested in the estimateCells function
  
  ###TEST1####
  nc1 <- ncol(getData(CIMseqSinglets_test, "counts"))
  nc2 <- ncol(getData(CIMseqMultiplets_test, "counts"))
  nc <- nc1 + nc2
  
  #run function
  output <- plotCountsData(CIMseqSinglets_test, CIMseqMultiplets_test)
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
  expect_true(ncol(output) == 7)
  expect_true(nrow(output) == nc)
})

##run test plotCountsERCC
test_that("check that plotCountsERCC outputs the expected result", {
  
  #run function
  output <- plotCountsERCC(CIMseqSinglets_test, CIMseqMultiplets_test)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotCountsMarkers
test_that("check that plotCountsMarkers outputs the expected result", {
    
  #run function
  markers <- c("RNR2", "NME7")
  output <- plotCountsMarkers(CIMseqSinglets_test, CIMseqMultiplets_test, markers)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
  expect_error(plotCountsMarkers(CIMseqSinglets_test, CIMseqMultiplets_test, "A"))
})

##run test plotUnsupervisedClass
test_that("check that plotUnsupervisedClass outputs the expected result", {
  
  #run function
  output <- plotUnsupervisedClass(CIMseqSinglets_test, CIMseqMultiplets_test)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotUnsupervisedMarkers
test_that("check that plotUnsupervisedMarkers outputs the expected result", {
  
  #run function
  output <- plotUnsupervisedMarkers(
    CIMseqSinglets_test, CIMseqMultiplets_test, markers = "ACTB"
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
  expect_error(plotUnsupervisedMarkers(
    CIMseqSinglets_test, CIMseqMultiplets_test, markers = NULL
  ))
  expect_error(plotUnsupervisedMarkers(
    CIMseqSinglets_test, CIMseqMultiplets_test, markers = "not a gene"
  ))
  
  #run function
  output <- plotUnsupervisedMarkers(
    CIMseqSinglets_test, CIMseqMultiplets_test,
    markers = c("ACTB", "GAPDH")
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

