context("All-classes")

test_that("check that accessors output the expected result", {
  expect_is(tsne(test_spUnsupervised), "matrix")
  expect_is(tsneMeans(test_spUnsupervised), "data.frame")
  expect_is(classification(test_spUnsupervised), "character")
  expect_is(uncertainty(test_spUnsupervised), "numeric")
  expect_is(selectInd(test_spUnsupervised), "integer")
})

test_that("check that replacement outputs the expected result", {
  expect_error(tsne(test_spUnsupervised) <- "a")
  expect_error(tsneMeans(test_spUnsupervised) <- "a")
  expect_error(classification(test_spUnsupervised) <- 1)
  expect_error(uncertainty(test_spUnsupervised) <- "a")
  expect_error(selectInd(test_spUnsupervised) <- "a")
  
  tsne(test_spUnsupervised) <- tsne(test_spUnsupervised)
  expect_is(test_spUnsupervised, "spUnsupervised")
  tsneMeans(test_spUnsupervised) <- tsneMeans(test_spUnsupervised)
  expect_is(test_spUnsupervised, "spUnsupervised")
  classification(test_spUnsupervised) <- classification(test_spUnsupervised)
  expect_is(test_spUnsupervised, "spUnsupervised")
  uncertainty(test_spUnsupervised) <- uncertainty(test_spUnsupervised)
  expect_is(test_spUnsupervised, "spUnsupervised")
  selectInd(test_spUnsupervised) <- selectInd(test_spUnsupervised)
  expect_is(test_spUnsupervised, "spUnsupervised")
})

