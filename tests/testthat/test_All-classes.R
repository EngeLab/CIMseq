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

test_that("check that CIMseqSinglets concatenation outputs the expected result", {
  #setup expected data
  counts <- cbind(
    getData(CIMseqSinglets_test, "counts"),
    getData(CIMseqSinglets_test, "counts")
  )
  counts.ercc <- cbind(
    getData(CIMseqSinglets_test, "counts.ercc"),
    getData(CIMseqSinglets_test, "counts.ercc")
  )
  dim.red <- rbind(
    getData(CIMseqSinglets_test, "dim.red"),
    getData(CIMseqSinglets_test, "dim.red")
  )
  classification <- c(
    getData(CIMseqSinglets_test, "classification"),
    getData(CIMseqSinglets_test, "classification")
  )
  
  expected <- CIMseqSinglets(counts, counts.ercc, dim.red, classification)
  
  #generate output
  output <- c(CIMseqSinglets_test, CIMseqSinglets_test)
  
  #test
  expect_identical(getData(expected, "counts"), getData(output, "counts"))
  expect_identical(getData(expected, "counts.log"), getData(output, "counts.log"))
  expect_identical(getData(expected, "counts.cpm"), getData(output, "counts.cpm"))
  expect_identical(getData(expected, "counts.ercc"), getData(output, "counts.ercc"))
  expect_identical(getData(expected, "dim.red"), getData(output, "dim.red"))
  expect_identical(getData(expected, "classification"), getData(output, "classification"))
  
  #TEST, 2 if genes differ
  CIMseqSinglets_test2 <- CIMseqSinglets(
    getData(CIMseqSinglets_test, "counts")[1:100, ],
    getData(CIMseqSinglets_test, "counts.ercc"),
    getData(CIMseqSinglets_test, "dim.red"),
    getData(CIMseqSinglets_test, "classification")
  )
  expect_error(c(CIMseqSinglets_test, CIMseqSinglets_test2))
  
  #TEST 3, ERCC differ
  CIMseqSinglets_test3 <- CIMseqSinglets(
    getData(CIMseqSinglets_test, "counts"),
    getData(CIMseqSinglets_test, "counts.ercc")[1:2, ],
    getData(CIMseqSinglets_test, "dim.red"),
    getData(CIMseqSinglets_test, "classification")
  )
  expect_error(c(CIMseqSinglets_test, CIMseqSinglets_test3))
  
  #TEST 4, different dim.red
  dr <- getData(CIMseqSinglets_test, "dim.red")
  colnames(dr) <- c("dim1", "dim2")
  CIMseqSinglets_test4 <- CIMseqSinglets(
    getData(CIMseqSinglets_test, "counts"),
    getData(CIMseqSinglets_test, "counts.ercc"),
    dr,
    getData(CIMseqSinglets_test, "classification")
  )
  expect_error(c(CIMseqSinglets_test, CIMseqSinglets_test4))
})

test_that("check that CIMseqMultiplets concatenation outputs the expected result", {
  #setup expected data
  counts <- cbind(
    getData(CIMseqMultiplets_test, "counts"),
    getData(CIMseqMultiplets_test, "counts")
  )
  counts.ercc <- cbind(
    getData(CIMseqMultiplets_test, "counts.ercc"),
    getData(CIMseqMultiplets_test, "counts.ercc")
  )
  features <- getData(CIMseqMultiplets_test, "features")
  expected <- CIMseqMultiplets(counts, counts.ercc, features)
  
  #generate output
  output <- c(CIMseqMultiplets_test, CIMseqMultiplets_test)
  
  #test
  expect_identical(getData(expected, "counts"), getData(output, "counts"))
  expect_identical(getData(expected, "counts.log"), getData(output, "counts.log"))
  expect_identical(getData(expected, "counts.cpm"), getData(output, "counts.cpm"))
  expect_identical(getData(expected, "counts.ercc"), getData(output, "counts.ercc"))
  expect_identical(getData(expected, "features"), getData(output, "features"))
  
  #TEST, 2 if genes differ
  CIMseqMultiplets_test2 <- CIMseqMultiplets(
    getData(CIMseqMultiplets_test, "counts")[1:100, ],
    getData(CIMseqMultiplets_test, "counts.ercc"),
    getData(CIMseqMultiplets_test, "features")
  )
  expect_error(c(CIMseqMultiplets_test, CIMseqMultiplets_test2))
  
  #TEST 3, ERCC differ
  CIMseqMultiplets_test3 <- CIMseqMultiplets(
    getData(CIMseqMultiplets_test, "counts"),
    getData(CIMseqMultiplets_test, "counts.ercc")[1:2, ],
    getData(CIMseqMultiplets_test, "features")
  )
  expect_error(c(CIMseqMultiplets_test, CIMseqMultiplets_test3))
  
  #TEST 4, different features
  CIMseqMultiplets_test4 <- CIMseqMultiplets(
    getData(CIMseqMultiplets_test, "counts"),
    getData(CIMseqMultiplets_test, "counts.ercc"),
    1:10
  )
  expect_warning(c(CIMseqMultiplets_test, CIMseqMultiplets_test4))
})

test_that("check that CIMseqSwarm concatenation outputs the expected result", {
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
  singletIdx <- c(getData(CIMseqSwarm_test, "singletIdx"))
  
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
  expect_identical(getData(expected, "fractions"), getData(output, "fractions"))
  expect_identical(getData(expected, "costs"), getData(output, "costs"))
  expect_identical(getData(expected, "convergence"), getData(output, "convergence"))
  expect_identical(getData(expected, "stats"), getData(output, "stats"))
  expect_identical(getData(expected, "singletIdx"), getData(output, "singletIdx"))
  expect_identical(getData(expected, "arguments"), getData(output, "arguments"))
  
  #TEST 2, non-matching classes
  f <- getData(CIMseqSwarm_test, "fractions")
  colnames(f) <- paste0("c", 1:ncol(f))
  CIMseqSwarm_test2 <- new(
    "CIMseqSwarm",
    fractions = f, 
    costs = getData(CIMseqSwarm_test, "costs"), 
    convergence = getData(CIMseqSwarm_test, "convergence"),
    stats = getData(CIMseqSwarm_test, "stats"),
    singletIdx = getData(CIMseqSwarm_test, "singletIdx"),
    arguments = getData(CIMseqSwarm_test, "arguments")
  )
  expect_error(c(CIMseqSwarm_test, CIMseqSwarm_test2))
  
  #TEST 3; non-matching singletIdx
  CIMseqSwarm_test3 <- new(
    "CIMseqSwarm",
    fractions = getData(CIMseqSwarm_test, "fractions"),
    costs = getData(CIMseqSwarm_test, "costs"), 
    convergence = getData(CIMseqSwarm_test, "convergence"),
    stats = getData(CIMseqSwarm_test, "stats"),
    singletIdx = getData(CIMseqSwarm_test, "singletIdx")[1:200],
    arguments = getData(CIMseqSwarm_test, "arguments")
  )
  expect_warning(c(CIMseqSwarm_test, CIMseqSwarm_test3))
})

test_that("check that CIMseqSwarm subsetting outputs the expected result", {
  #setup expected data
  s <- c("m.NJB00204.G04", "m.NJB00204.D07")
  fractions <- getData(CIMseqSwarm_test, "fractions")[s, ]
  costs <- getData(CIMseqSwarm_test, "costs")[s]
  convergence <- getData(CIMseqSwarm_test, "convergence")[s]
  stats <- filter(getData(CIMseqSwarm_test, "stats"), sample %in% s)
  singletIdx <- getData(CIMseqSwarm_test, "singletIdx")
  arguments <- getData(CIMseqSwarm_test, "arguments")
  
  expected <- new(
    "CIMseqSwarm",
    fractions = fractions, costs = costs, convergence = convergence,
    stats = stats, singletIdx = singletIdx, arguments = arguments
  )
  
  #generate output
  output <- filterSwarm(CIMseqSwarm_test, s)
  
  #test
  expect_identical(getData(expected, "fractions"), getData(output, "fractions"))
  expect_identical(getData(expected, "costs"), getData(output, "costs"))
  expect_identical(getData(expected, "convergence"), getData(output, "convergence"))
  expect_identical(getData(expected, "stats"), getData(output, "stats"))
  expect_identical(getData(expected, "singletIdx"), getData(output, "singletIdx"))
  expect_identical(getData(expected, "arguments"), getData(output, "arguments"))
})

