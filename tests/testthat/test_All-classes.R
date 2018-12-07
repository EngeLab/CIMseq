context("All-classes")

test_that("check that accessors output the expected result", {
  expect_is(getData(CIMseqSinglets_test, "counts"), "matrix")
  expect_is(getData(CIMseqSinglets_test, "counts.log"), "matrix")
  expect_is(getData(CIMseqSinglets_test, "counts.cpm"), "matrix")
  expect_is(getData(CIMseqSinglets_test, "counts.ercc"), "matrix")
  expect_is(getData(CIMseqSinglets_test, "dim.red"), "matrix")
  expect_is(getData(CIMseqSinglets_test, "classification"), "character")
  
  expect_is(getData(CIMseqMultiplets_test, "counts"), "matrix")
  expect_is(getData(CIMseqMultiplets_test, "counts.log"), "matrix")
  expect_is(getData(CIMseqMultiplets_test, "counts.cpm"), "matrix")
  expect_is(getData(CIMseqMultiplets_test, "counts.ercc"), "matrix")
  expect_is(getData(CIMseqMultiplets_test, "features"), "integer")
  
  expect_is(getData(CIMseqSwarm_test, "fractions"), "matrix")
  expect_is(getData(CIMseqSwarm_test, "costs"), "numeric")
  expect_is(getData(CIMseqSwarm_test, "convergence"), "character")
  expect_is(getData(CIMseqSwarm_test, "stats"), "tbl_df")
  expect_is(getData(CIMseqSwarm_test, "singletIdx"), "list")
  expect_is(getData(CIMseqSwarm_test, "arguments"), "tbl_df")
})

test_that("check that replacement outputs the expected result", {
  expect_error(getData(CIMseqSinglets_test, "counts") <- "a")
  expect_error(getData(CIMseqSinglets_test, "counts.log") <- "a")
  expect_error(getData(CIMseqSinglets_test, "counts.cpm") <- "a")
  expect_error(getData(CIMseqSinglets_test, "counts.ercc") <- "a")
  expect_error(getData(CIMseqSinglets_test, "dim.red") <- "a")
  expect_error(getData(CIMseqSinglets_test, "classification") <- 1)
  
  expect_error(getData(CIMseqMultiplets_test, "counts") <- "a")
  expect_error(getData(CIMseqMultiplets_test, "counts.log") <- "a")
  expect_error(getData(CIMseqMultiplets_test, "counts.cpm") <- "a")
  expect_error(getData(CIMseqMultiplets_test, "counts.ercc") <- "a")
  expect_error(getData(CIMseqMultiplets_test, "features") <- "a")
  
  tmp <- CIMseqSinglets_test
  getData(tmp, "counts") <- getData(CIMseqSinglets_test, "counts")
  expect_is(tmp, "CIMseqSinglets")

  tmp <- CIMseqSinglets_test
  getData(tmp, "counts.log") <- function(){}
  expect_is(tmp, "CIMseqSinglets")

  tmp <- CIMseqSinglets_test
  getData(tmp, "counts.cpm") <- function(){}
  expect_is(tmp, "CIMseqSinglets")

  tmp <- CIMseqSinglets_test
  getData(tmp, "counts.ercc") <- getData(CIMseqSinglets_test, "counts.ercc")
  expect_is(tmp, "CIMseqSinglets")

  tmp <- CIMseqSinglets_test
  getData(tmp, "dim.red") <- getData(CIMseqSinglets_test, "dim.red")
  expect_is(tmp, "CIMseqSinglets")

  tmp <- CIMseqSinglets_test
  getData(tmp, "classification") <- getData(CIMseqSinglets_test, "classification")
  expect_is(tmp, "CIMseqSinglets")

  tmp <- CIMseqMultiplets_test
  getData(tmp, "counts") <- getData(CIMseqMultiplets_test, "counts")
  expect_is(tmp, "CIMseqMultiplets")

  tmp <- CIMseqMultiplets_test
  getData(tmp, "counts.log") <- function(){}
  expect_is(tmp, "CIMseqMultiplets")

  tmp <- CIMseqMultiplets_test
  getData(tmp, "counts.cpm") <- function(){}
  expect_is(tmp, "CIMseqMultiplets")

  tmp <- CIMseqMultiplets_test
  getData(tmp, "counts.ercc") <- getData(CIMseqMultiplets_test, "counts.ercc")
  expect_is(tmp, "CIMseqMultiplets")

  tmp <- CIMseqMultiplets_test
  getData(tmp, "features") <- getData(CIMseqMultiplets_test, "features")
  expect_is(tmp, "CIMseqMultiplets")
})

test_that("check that concatenation outputs the expected result", {
  #setup expected data
  fractions <- rbind(
    getData(CIMseqSwarm_test, "fractions"),
    getData(CIMseqSwarm_test, "fractions")
  )
  costs <- c(
    getData(CIMseqSwarm_test, "costs"), 
    getData(CIMseqSwarm_test, "costs")
  )
  convergence <- c(
    getData(CIMseqSwarm_test, "convergence"), 
    getData(CIMseqSwarm_test, "convergence")
  )
  stats <- rbind(
    getData(CIMseqSwarm_test, "stats"),
    getData(CIMseqSwarm_test, "stats")
  )
  singletIdx <- c(
    getData(CIMseqSwarm_test, "singletIdx"), 
    getData(CIMseqSwarm_test, "singletIdx")
  )
  arguments <- rbind(
    getData(CIMseqSwarm_test, "arguments"),
    getData(CIMseqSwarm_test, "arguments")
  )
  
  expected <- new("CIMseqSwarm",
    fractions = fractions, costs = costs, convergence = convergence,
    stats = stats, singletIdx = singletIdx, arguments = arguments
  )
  
  #generate output
  output <- c(CIMseqSwarm_test, CIMseqSwarm_test)
  
  #test
  expect_identical(expected, output)
})

