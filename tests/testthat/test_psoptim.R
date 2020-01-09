context("psoptim")


test_that("check that swarmInit runs without issue", {
  output <- swarmInit(CIMseqSinglets_test, 2)
  expect_true(nrow(output) == 3)
  expect_true(ncol(output) == 6)
  expect_is(output, "matrix")
  expect_true(mode(output) == "numeric")
  expect_error(swarmInit(CIMseqSinglets_test, 1))
})

