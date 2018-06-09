context("spSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])
uObj <- testUns
sObj <- testSwa


#Function to check if all elements in a vector are identical
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

##run test getMultipletsForEdge
test_that("check that getMultipletsForEdge outputs the expected result", {
  
    ###TEST1####
    #setup expected data
    #A1 and B1 should have an edge
    #I1 and J1 should have an edge
    expected1 <- tibble::tibble(
        multiplet = "m.A1B1.341",
        from = "A1",
        to = "B1"
    )
    expected2 <- tibble::tibble(
        multiplet = "m.C1D1.342",
        from = "C1",
        to = "D1"
    )
    expected3 <- tibble::tibble(
        multiplet = c("m.A1B1.341", "m.C1D1.342"),
        from = c("A1", "C1"),
        to = c("B1", "D1")
    )
    
    #run function
    output1 <- getMultipletsForEdge(sObj, 1/4, data.frame("A1", "B1"))
    output2 <- getMultipletsForEdge(sObj, 1/4, data.frame("C1", "D1"))
    output3 <- getMultipletsForEdge(
      sObj, 1/4, data.frame(c("A1", "C1"), c("B1", "D1"))
    )
    
    #test
    expect_identical(output1, expected1)
    expect_identical(output2, expected2)
    expect_identical(output3, expected3)
    
    ###TEST2####
    #prepare normal input data
    tmp <- sObj
    tmp@spSwarm <-
        data.frame(
            A1 = c(0.3, 0.3),
            B1 = c(0.3, 0.3),
            C1 = c(0.3, 0),
            D1 = c(0, 0.3),
            row.names=c("m.A1B1C1a", "m.A1B1C1b")
        )
    
    #setup expected data
    expected <- tibble::tibble(
        multiplet = rep(c("m.A1B1C1a", "m.A1B1C1b"), 3),
        from = c(rep("A1", 4), rep("B1", 2)),
        to = c(rep("B1", 2), rep(c("C1", "D1"), 2))
    )
    
    
    #run function
    output <- getMultipletsForEdge(
        tmp,
        1/4,
        data.frame(
            c("A1", "A1", "A1", "B1", "B1", "C1"),
            c("B1", "C1", "D1", "C1", "D1", "D1")
        )
    )
    
    #test
    expect_identical(output, expected)
})

##run test getEdgesForMultiplet
test_that("check that getEdgesForMultiplet outputs the expected result", {
  
    ###TEST1####
    #setup expected data
    expected1 <- tibble::tibble(
        multiplet = "m.A1B1",
        from = "A1",
        to = "B1"
    )
    expected2 <- tibble::tibble(
        multiplet = "m.C1D1",
        from = "C1",
        to = "D1"
    )

    #run function
    output1 <- getEdgesForMultiplet(sObj, 1/4, 'm.A1B1')
    output2 <- getEdgesForMultiplet(sObj, 1/4, 'm.C1D1')

    #test
    expect_identical(output1, expected1)
    expect_identical(output2, expected2)
})

################################################################################
#                                                                              #
#                                C++ functions                                 #
#                                                                              #
################################################################################

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
  output <- sampleSinglets(classes)
  
  #test
  expect_true(length(output) == 3)
  expect_true(output[1] %in% a)
  expect_true(output[2] %in% b)
  expect_true(output[3] %in% c)
})

##run test subsetSingletsEigen
test_that("check that subsetSingletsEigen outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 10), ncol = 10)
  idx <- c(0, 7)
  
  #setup expected data
  a <- 1
  b <- 8
  
  #run function
  output <- subsetSingletsEigen(singlets, idx)
  
  #test
  expect_true(all(output[, 1] == a))
  expect_true(all(output[, 2] == b))
})

##run test adjustAccordingToFractionsEigen
test_that("check that adjustAccordingToFractionsEigen outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 10), ncol = 10) * 0.5
  fractions <- rep(c(0.1, 0.9), each = 5)
  
  #setup expected data
  expected <- t(t(singlets) * fractions)
  
  #run function
  output <- adjustAccordingToFractionsEigen(fractions, singlets)
  
  #test
  expect_identical(output, expected)
})

##run test multipletSumsEigen
test_that("check that multipletSumsEigen outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(c(rep(0.5, 10), rep(2:10, each = 10)), ncol = 10)
  
  #setup expected data
  ##note that R rounds down with 0.5 i.e. round(54.5) gives 54 whereas Armadillo
  ##gives 55. Due to this, ceiling is used when generating expected.
  expected <- ceiling(rowSums(singlets))
  
  #run function
  output <- multipletSumsEigen(singlets)
  
  #test
  expect_identical(output, expected)
})

##run test poissonSample
test_that("check that poissonSample outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 1000), ncol = 1000)
  rs <- rowSums(singlets)
  
  #setup expected data
  expected <- rpois(length(rs), rs)
  
  #run function
  output <- poissonSample(matrix(rs, ncol = 1))
  
  #test
  expect_true(all.equal(mean(output), mean(expected), tolerance = 100))
})

##run test cpmC
test_that("check that cpmC outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  counts <- matrix(sample(seq(0.01, 0.001, -0.00001), 20) + 1, ncol = 2)
  
  #setup expected data
  expected <- t(t(counts) / colSums(counts) * 10^6 + 1)
  
  #run function
  output <- cpmC(counts)
  
  #test
  expect_equal(output, expected)
  expect_error(cpmC(matrix(1:10, ncol = 2)))
})

##run test generateSyntheticMultipletsEigen
test_that("check that generateSyntheticMultipletsEigen outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 10), ncol = 10)
  fractions <- c(0.5, 1)
  classes <- rep(LETTERS[1:2], each = 5)
  
  #run function
  output <- generateSyntheticMultipletsEigen(singlets, classes, fractions, 100)
  
  #test
  expect_silent(generateSyntheticMultipletsEigen(singlets, classes, fractions, 1))
  expect_silent(generateSyntheticMultipletsEigen(singlets, classes, fractions, 2))
  expect_false(zero_range(output[1, ]))
})

##run test calculateCostDensity
test_that("check that calculateCostDensity outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  set.seed(3289)
  oneMultiplet <- as.integer(runif(10, 1, 10))
  syntheticMultiplets <- matrix(runif(20, 1, 10), ncol = 2)
  
  #setup expected data
  expected <- dpois(oneMultiplet, syntheticMultiplets)
  
  #run function
  output <- calculateCostDensity(oneMultiplet, syntheticMultiplets)
  
  #test
  expect_identical(output, expected)
})

##run test calculateLogRowMeans
test_that("check that calculateLogRowMeans outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  densities <- matrix(sample(seq(0.001, 0.1, 0.0001), 20), ncol = 2)
  
  #setup expected data
  expected <- log10(rowMeans(densities))
  
  #run function
  output <- calculateLogRowMeans(densities)
  
  #test
  expect_identical(output, expected)
})

##run test fixNegInf
test_that("check that fixNegInf outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  values <- c(-Inf, rnorm(20), -Inf)
  
  #setup expected data
  tmp <- values
  tmp[is.infinite(tmp)] <- -323.0052
  expected <- tmp
  
  #run function
  output <- fixNegInf(values)
  
  #test
  expect_identical(output, expected)
})

##run test costNegSum
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
