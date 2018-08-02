context("spSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])
uObj <- testUns
sObj <- testSwa


#Function to check if all elements in a vector are identical
has_zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

##run test getMultipletsForEdge
# test_that("check that getMultipletsForEdge outputs the expected result", {
#   
#     ###TEST1####
#     #setup expected data
#     #A1 and B1 should have an edge
#     #I1 and J1 should have an edge
#     expected1 <- tibble::tibble(
#         multiplet = "m.A1B1.341",
#         from = "A1",
#         to = "B1"
#     )
#     expected2 <- tibble::tibble(
#         multiplet = "m.C1D1.342",
#         from = "C1",
#         to = "D1"
#     )
#     expected3 <- tibble::tibble(
#         multiplet = c("m.A1B1.341", "m.C1D1.342"),
#         from = c("A1", "C1"),
#         to = c("B1", "D1")
#     )
#     
#     #run function
#     output1 <- getMultipletsForEdge(sObj, 0, data.frame("A1", "B1"))
#     output2 <- getMultipletsForEdge(sObj, 1/4, data.frame("C1", "D1"))
#     output3 <- getMultipletsForEdge(
#       sObj, 1/4, data.frame(c("A1", "C1"), c("B1", "D1"))
#     )
#     
#     #test
#     expect_identical(output1, expected1)
#     expect_identical(output2, expected2)
#     expect_identical(output3, expected3)
#     
#     ###TEST2####
#     #prepare normal input data
#     tmp <- sObj
#     tmp@spSwarm <-
#         data.frame(
#             A1 = c(0.3, 0.3),
#             B1 = c(0.3, 0.3),
#             C1 = c(0.3, 0),
#             D1 = c(0, 0.3),
#             row.names=c("m.A1B1C1a", "m.A1B1C1b")
#         )
#     
#     #setup expected data
#     expected <- tibble::tibble(
#         multiplet = rep(c("m.A1B1C1a", "m.A1B1C1b"), 3),
#         from = c(rep("A1", 4), rep("B1", 2)),
#         to = c(rep("B1", 2), rep(c("C1", "D1"), 2))
#     )
#     
#     
#     #run function
#     output <- getMultipletsForEdge(
#         tmp,
#         1/4,
#         data.frame(
#             c("A1", "A1", "A1", "B1", "B1", "C1"),
#             c("B1", "C1", "D1", "C1", "D1", "D1")
#         )
#     )
#     
#     #test
#     expect_identical(output, expected)
# })
# 
# ##run test getEdgesForMultiplet
# test_that("check that getEdgesForMultiplet outputs the expected result", {
#   
#     ###TEST1####
#     #setup expected data
#     expected1 <- tibble::tibble(
#         multiplet = "m.A1B1",
#         from = "A1",
#         to = "B1"
#     )
#     expected2 <- tibble::tibble(
#         multiplet = "m.C1D1",
#         from = "C1",
#         to = "D1"
#     )
# 
#     #run function
#     output1 <- getEdgesForMultiplet(sObj, 1/4, 'm.A1B1')
#     output2 <- getEdgesForMultiplet(sObj, 1/4, 'm.C1D1')
# 
#     #test
#     expect_identical(output1, expected1)
#     expect_identical(output2, expected2)
# })

################################################################################
#                                                                              #
#                                C++ functions                                 #
#                                                                              #
################################################################################

context("generateSyntheticMultiplets")

################################################################################
##          ARMADILLO FUNCTIONS TO GENERATE SYNTHETIC MULTIPLETS              ##
################################################################################

##run test normalizeFractions
test_that("check that normalizeFractions outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  fractions <- c(0.1, 0.2)
  
  #setup expected data
  expected <- matrix(fractions / sum(fractions))
  
  #run function
  output <- normalizeFractions(fractions)
  
  #test
  expect_equal(output, expected)
})

##run test sampleSinglets
test_that("check that sampleSinglets outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  classes <- rep(LETTERS[1:3], each = 2)
  
  #setup expected data
  a <- c(0, 1)
  b <- c(2, 3)
  c <- c(4, 5)
  
  #run function
  output <- sampleSinglets(classes) #returns a matrix
  
  #test
  expect_true(length(output) == 3)
  expect_true(output[1] %in% a)
  expect_true(output[2] %in% b)
  expect_true(output[3] %in% c)
})

##run test subsetSinglets
test_that("check that subsetSinglets outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 10), ncol = 10)
  idx <- c(0, 7)
  
  #setup expected data
  a <- 1
  b <- 8
  
  #run function
  output <- lapply(1:10, function(x) {
    set.seed(30932 + x)
    subsetSinglets(singlets, idx)
  }) %>%
    do.call("rbind", .)
  
  #test
  expect_true(all(output[, 1] == a))
  expect_true(all(output[, 2] == b))
  
  ###TEST2####
  #prepare normal input data
  classes <- rep(letters[1:4], each = 5)
  a <- rep(1, 10)
  b <- rep(2, 10)
  c <- rep(3, 10)
  d <- rep(4, 10)
  singlets <- cbind(
    a, a, a, a, a, b, b, b, b, b, c, c, c, c, c, d, d, d, d, d
  )
  rownames(singlets) <- LETTERS[1:nrow(singlets)]
  n <- 10
  
  #setup expected data
  expected <- matrix(rep(1:4, each = 10), nrow = 10)
  storage.mode(expected) <- "numeric"
  
  #run function
  set.seed(82390)
  idx <- purrr::map(1:n, ~sampleSinglets(classes))
  
  #check the indexes first
  idxCheck <- purrr::map_dfc(idx, function(i) tibble::as_tibble(classes[i + 1]))
  expect_true(all("a" == unname(unlist(idxCheck[1, ]))))
  expect_true(all("b" == unname(unlist(idxCheck[2, ]))))
  expect_true(all("c" == unname(unlist(idxCheck[3, ]))))
  expect_true(all("d" == unname(unlist(idxCheck[4, ]))))
  
  #run subsetSinglets
  output <- purrr::map(idx, ~subsetSinglets(singlets, .x))
  
  #test
  expect_true(all(purrr::map_lgl(output, ~identical(.x, expected))))
})

##run test adjustAccordingToFractions
test_that("check that adjustAccordingToFractions outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 10), ncol = 10) * 0.5
  fractions <- rep(c(0.1, 0.9), each = 5)
  
  #setup expected data
  expected <- t(t(singlets) * fractions)
  
  #run function
  output <- adjustAccordingToFractions(fractions, singlets)
  
  #test
  expect_identical(output, expected)
})

##run test multipletSums
test_that("check that multipletSums outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(c(rep(0.5, 10), rep(2:10, each = 10)), ncol = 10)
  
  #setup expected data
  ##note that R rounds down with 0.5 i.e. round(54.5) gives 54 whereas Armadillo
  ##gives 55. Due to this, ceiling is used when generating expected.
  expected <- matrix(ceiling(rowSums(singlets)))
  
  #run function
  output <- multipletSums(singlets)
  
  #test
  expect_identical(output, expected)
})

##run test vecToMat
test_that("check that vecToMat outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  vec <- c(
    rep(1, 10),
    rep(2, 10),
    rep(3, 10)
  )
  names(vec) <- c(
    rep("a", 10),
    rep("b", 10),
    rep("c", 10)
  )
  #setup expected data
  expected <- matrix(vec, nrow = 3, byrow = TRUE)
  
  #run function
  output <- vecToMat(vec, 3, 10)
  
  #test
  expect_identical(output, expected)
})


################################################################################
##                        FUNCTIONS TO CALCULATE COST ARMA                    ##
################################################################################

##run test calculateCostDensityArma
test_that("check that calculateCostDensity outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  set.seed(3289)
  oneMultiplet <- as.integer(runif(10, 1, 10))
  syntheticMultiplets <- matrix(runif(20, 1, 10) + 0.1, ncol = 2)
  
  #setup expected data
  expected <- dpois(oneMultiplet, syntheticMultiplets)
  
  #run function
  storage.mode(oneMultiplet) <- "numeric"
  output <- calculateCostDensity(oneMultiplet, syntheticMultiplets)
  
  #test
  expect_identical(output, expected)
})

##run test calculateLogRowMeansArma
test_that("check that calculateLogRowMeans outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  densities <- matrix(sample(seq(0.001, 0.1, 0.0001), 20), ncol = 2)
  
  #setup expected data
  expected <- matrix(log10(rowMeans(densities)))
  
  #run function
  output <- calculateLogRowMeans(densities)
  
  #test
  expect_identical(output, expected)
})

##run test fixNegInfArma
test_that("check that fixNegInf outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  values <- c(-Inf, rnorm(20), -Inf)
  
  #setup expected data
  tmp <- values
  tmp[is.infinite(tmp)] <- -323.0052
  expected <- matrix(tmp)
  
  #run function
  output <- fixNegInf(values)
  
  #test
  expect_identical(output, expected)
})

##run test costNegSumArma
test_that("check that costNegSum outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  values <- rnorm(20)
  
  #setup expected data
  expected <- -sum(values)
  
  #run function
  output <- costNegSum(values)
  
  #test
  expect_true(all.equal(output, expected))
})

##run test .subsetSinglets
test_that("check that .subsetSinglets outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  classes <- rep(letters[1:4], each = 5)
  a <- rep(1, 10)
  b <- rep(2, 10)
  c <- rep(3, 10)
  d <- rep(4, 10)
  singlets <- cbind(a, a, a, a, a, b, b, b, b, b, c, c, c, c, c, d, d, d, d, d)
  rownames(singlets) <- LETTERS[1:nrow(singlets)]
  n <- 10
  
  #setup expected data
  nc <- length(unique(classes))
  nr <- nrow(singlets) * n
  rn <- rep(LETTERS[1:nrow(singlets)], each = n)
  col1 <- rep(1, n * nrow(singlets))
  col2 <- rep(2, n * nrow(singlets))
  col3 <- rep(3, n * nrow(singlets))
  col4 <- rep(4, n * nrow(singlets))
  
  #run function
  set.seed(82390)
  output <- .subsetSinglets(classes, singlets, n)
  
  #test
  expect_equal(nr, nrow(output))
  expect_equal(nc, ncol(output))
  expect_equal(col1, unname(output[, 1]))
  expect_equal(col2, unname(output[, 2]))
  expect_equal(col3, unname(output[, 3]))
  expect_equal(col4, unname(output[, 4]))
  expect_identical(rn, rownames(output))
})

test_that("check that calculateCost and cost give identical results", {
  
  ###TEST1####
  #prepare normal input data
  singletSubset <- getData(sObj, "syntheticMultiplets")
  m <- "m.A1B1.341"
  oneMultiplet <- ceiling(getData(cObjMul, "counts.cpm")[, m])
  fractions <- as.numeric(getData(sObj, "spSwarm")[m, ])
  n <- getData(sObj, "arguments")$nSyntheticMultiplets
  
  #expected
  expected <- cost(oneMultiplet, singletSubset, fractions, n)
  
  #output
  output <- calculateCost(oneMultiplet, singletSubset, fractions, n)
  
  #test
  expect_equal(expected, output)
})

