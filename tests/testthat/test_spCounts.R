context("spCounts")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[, s], testErcc[, s])
cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])

test_that("check that the .norm.counts function outputs the expected result", {
    
  ###TEST1####
  #prepare normal input data
  input <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  
  #setup expected data
  expected <- matrix(
    c(0, 1000000, 1000000, 0),
    nrow = 2,
    ncol = 2
  )
  
  #run function
  output <- .norm.counts(input)
  
  #test
  expect_true(all.equal(expected, output))
})

test_that("check that the .inputCheckCounts outputs the expected result", {
    
    ###TEST1####
    #non-conforming counts and counts.ercc
    #setup input
    counts <- matrix(1:10, ncol = 10)
    counts.ercc <- matrix(1:11, ncol = 11)
    
    #test
    expect_message(.inputCheckCounts(counts, counts.ercc))
    
    ###TEST2####
    #NA's present
    #setup input
    counts <- matrix(NA, ncol = 10)
    counts.ercc <- matrix(NA, ncol = 11)
    
    #test
    expect_message(.inputCheckCounts(counts, counts.ercc))
})

test_that("check that estimateCells outputs the expected result", {
    
  ###TEST1####
  #setup expected data
  expected1 <- tibble::tibble(
    sampleType=c(rep("Singlet", 340), rep("Multiplet", 2))
  )
  
  expected2 <- tibble::tibble(
    frac.ercc = 1,
    cellNumberMin = 1,
    cellNumberMedian = 1,
    cellNumberMax = 2
  )
  
  expected3 <- tibble::tibble(
    frac.ercc = 0,
    cellNumberMin = 2,
    cellNumberMedian = 3,
    cellNumberMax = 3
  )
  
  #run function
  output <- estimateCells(cObjSng, cObjMul)
  firstRow <- output %>%
    dplyr::slice(1) %>%
    dplyr::select(frac.ercc:cellNumberMax) %>%
    round()
  lastRow <- output %>%
    dplyr::slice(nrow(output)) %>%
    dplyr::select(frac.ercc:cellNumberMax) %>%
    round()
  
  #test
  expect_identical(expected1, dplyr::select(output, `sampleType`))
  expect_identical(expected2, firstRow)
  expect_identical(expected3, lastRow)
  expect_type(output$sampleName, "character")
  expect_type(output$sampleType, "character")
  expect_type(output$frac.ercc, "double")
  expect_type(output$cellNumberMin, "double")
  expect_type(output$cellNumberMedian, "double")
  expect_type(output$cellNumberMax, "double")
})

test_that("check that estimateCells gives error with all 0 ercc", {
  ercc <- getData(cObjSng, "counts.ercc")
  ercc[ ,1] <- rep(0, nrow(ercc))
  cObjSng@counts.ercc <- ercc
  expect_warning(estimateCells(cObjSng, cObjMul))
})
