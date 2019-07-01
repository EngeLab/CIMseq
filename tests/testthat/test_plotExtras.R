context("plotExtras")

##run test coloursFromTargets
test_that("check that coloursFromTargets outputs the expected result", {
  
  ###TEST1####
  #run function
  singlets <- getData(CIMseqSinglets_test, "counts.cpm")
  output <- coloursFromTargets(col40(), singlets, rownames(singlets)[1:3])
  
  #test
  expect_is(output, c("tbl_df", "tbl", "data.frame"))
  expect_type(output, "list")
  
  ###TEST2####
  counts <- getData(CIMseqSinglets_test, "counts")
  output <- coloursFromTargets(letters[1:2], counts, "A")
  expect_identical(tibble(Sample = colnames(counts)), output)
  
  ###TEST3####
  counts[1, ] <- 0
  expect_error(coloursFromTargets(letters[1:2], counts, c("RNR2", "ACTB")))
})

test_that("check that plotData outputs the expected result", {
  plot <- ggplot(mtcars) + geom_point(aes(mpg, drat))
  expect_identical(as_tibble(mtcars), plotData(plot))
})

test_that("check that convertToERCC outputs the expected result", {
  expected <- c(
    1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 
    0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 
    1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 2, 1, 1, 
    1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0
  )
  ercc <- estimateCells(CIMseqSinglets_test, CIMseqMultiplets_test)
  output <- convertToERCC(
    ercc$estimatedCellNumber, CIMseqSinglets_test, CIMseqMultiplets_test
  )
  expect_identical(expected, round(output))
})