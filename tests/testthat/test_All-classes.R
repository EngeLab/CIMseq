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
  expect_is(getData(CIMseqSwarm_test, "arguments"), "list")
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
  
  getData(CIMseqSinglets_test, "counts") <- getData(CIMseqSinglets_test, "counts")
  expect_is(CIMseqSinglets_test, "CIMseqSinglets")
  getData(CIMseqSinglets_test, "counts.log") <- function(){}
  expect_is(CIMseqSinglets_test, "CIMseqSinglets")
  getData(CIMseqSinglets_test, "counts.cpm") <- function(){}
  expect_is(CIMseqSinglets_test, "CIMseqSinglets")
  getData(CIMseqSinglets_test, "counts.ercc") <- getData(CIMseqSinglets_test, "counts.ercc")
  expect_is(CIMseqSinglets_test, "CIMseqSinglets")
  getData(CIMseqSinglets_test, "dim.red") <- getData(CIMseqSinglets_test, "dim.red")
  expect_is(CIMseqSinglets_test, "CIMseqSinglets")
  getData(CIMseqSinglets_test, "classification") <- getData(CIMseqSinglets_test, "classification")
  expect_is(CIMseqSinglets_test, "CIMseqSinglets")
  
  getData(CIMseqMultiplets_test, "counts") <- getData(CIMseqMultiplets_test, "counts")
  expect_is(CIMseqMultiplets_test, "CIMseqMultiplets")
  getData(CIMseqMultiplets_test, "counts.log") <- function(){}
  expect_is(CIMseqMultiplets_test, "CIMseqMultiplets")
  getData(CIMseqMultiplets_test, "counts.cpm") <- function(){}
  expect_is(CIMseqMultiplets_test, "CIMseqMultiplets")
  getData(CIMseqMultiplets_test, "counts.ercc") <- getData(CIMseqMultiplets_test, "counts.ercc")
  expect_is(CIMseqMultiplets_test, "CIMseqMultiplets")
  getData(CIMseqMultiplets_test, "features") <- getData(CIMseqMultiplets_test, "features")
  expect_is(CIMseqMultiplets_test, "CIMseqMultiplets")
})

