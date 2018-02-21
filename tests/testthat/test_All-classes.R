#context("All-classes")

uObj <- testUns

test_that("check that accessors output the expected result", {
  expect_is(tsne(uObj), "matrix")
  expect_is(tsneMeans(uObj), "data.frame")
  expect_is(groupMeans(uObj), "matrix")
  expect_is(classification(uObj), "character")
  expect_is(uncertainty(uObj), "numeric")
  expect_is(selectInd(uObj), "integer")
})

test_that("check that replacement outputs the expected result", {
  expect_error(tsne(uObj) <- "a")
  expect_error(tsneMeans(uObj) <- "a")
  expect_error(groupMeans(uObj) <- "a")
  expect_error(classification(uObj) <- 1)
  expect_error(uncertainty(uObj) <- "a")
  expect_error(selectInd(uObj) <- "a")
})

