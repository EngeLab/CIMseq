context("tidyFunctions")

test_that("check that namedListToTibble outputs the expected result", {
  expected <- tibble::tibble(names = rep(letters[1:2], each = 10), variables = 1:20)
  output <- namedListToTibble(list(a = 1:10, b = 11:20))
  expect_identical(output, expected)
  
  input <- list(
    list(a = 1:10, b = 11:20),
    list(c = 21:30, d = 31:40)
  )
  expect_error(namedListToTibble(input))
  
  vec <- 1:3
  names(vec) <- letters[1:3]
  output <- namedListToTibble(list(a = vec, b = vec))
  expected <- tibble(
    names = rep(letters[1:2], each = 3),
    inner.names = rep(letters[1:3], 2),
    variables = rep(1:3, 2)
  )
  expect_identical(expected, output)
})

test_that("check that matrix_to_tibble outputs the expected result", {
  input <- matrix(1:10, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
  expected <- tibble::tibble(rowname = letters[1:5], A = 1:5, B = 6:10)
  output <- matrix_to_tibble(input)
  expect_identical(output, expected)
  
  input <- matrix(1:10, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
  expected <- tibble::tibble(A = 1:5, B = 6:10)
  output <- matrix_to_tibble(input, drop = TRUE)
  expect_identical(output, expected)
  
  input <- data.frame(A = 1:5, B = 6:10, row.names = letters[1:5])
  expect_error(matrix_to_tibble(input))
  
  input <- matrix(1:10, ncol = 2)
  output <- matrix_to_tibble(input)
  expected <- tibble(rowname = 1:5, V1 = 1:5, V2 = 6:10)
  expect_identical(expected, output)
})

test_that("check that tidySwarm outputs the expected result", {
  expected <- tibble(
    sample = c("m.NJB00204.G04", "m.NJB00204.D07", "m.NJB00204.A02"),
    costs = round(c(5594.8879991901, 8876.00099961946, 5532.41954319017)),
    convergence = rep("Maximal number of iterations reached.", 3),
    A375 = round(c(1, 0, 1)),
    HCT116 = round(c(0, 0, 0)),
    HOS = round(c(0, 1, 0))
  )
  output <- tidySwarm(CIMseqSwarm_test) %>%
    unnest(cols = .data$fractions) %>%
    mutate_if(is.numeric, round)
  expect_identical(expected, output)
})