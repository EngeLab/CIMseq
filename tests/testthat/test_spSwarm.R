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
    output1 <- getMultipletsForEdge(sObj, 0, data.frame("A1", "B1"))
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

context("generateSyntheticMultiplets")

################################################################################
##          ARMADILLO FUNCTIONS TO GENERATE SYNTHETIC MULTIPLETS              ##
################################################################################

##run test normalizeFractionsArma
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

##run test sampleSingletsArma
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

##run test subsetSingletsArma
test_that("check that subsetSinglets outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- matrix(rep(1:10, each = 10), ncol = 10)
  idx <- c(0, 7)
  
  #setup expected data
  a <- 1
  b <- 8
  
  #run function
  output <- subsetSinglets(singlets, idx)
  
  #test
  expect_true(all(output[, 1] == a))
  expect_true(all(output[, 2] == b))
})

##run test adjustAccordingToFractionsArma
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

##run test multipletSumsArma
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
  singlets <- cbind(a, a, a, a, a, b, b, b, b, b, c, c, c, c, c, c, d, d, d, d, d)
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
  output <- .subsetSinglets(classes, singlets, n)
  
  #test
  expect_equal(nr, nrow(output))
  expect_equal(nc, ncol(output))
  expect_equal(col1, unname(output[, 1]))
  expect_equal(col2, unname(output[, 2]))
  expect_equal(col3, unname(output[, 3]))
  expect_equal(col4, unname(output[, 4]))
  expect_identical(rn, rownames(output))
  
  ##TEST2
  classes <- c("SI.Stem", "C.Stem", "Blood", "C.Stem", "C.Colonocyte", "SI.Enterocyte", 
               "C.Goblet", "SI.Goblet", "C.Goblet", "C.Stem", "SI.Enterocyte", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", "C.Goblet", 
               "SI.Stem", "C.Colonocyte", "C.Goblet", "C.Stem", "C.Stem", "SI.Enterocyte", 
               "SI.Enterocyte", "C.Goblet", "C.Goblet", "C.Stem", "SI.Enterocyte", 
               "SI.Stem", "C.Colonocyte", "C.Colonocyte", "C.Goblet", "C.Stem", 
               "SI.Stem", "SI.Tufft", "C.Stem", "C.Goblet", "C.Stem", "C.Goblet", 
               "C.Stem", "SI.Stem", "C.Colonocyte", "SI.Enterocyte", "C.Stem", 
               "Blood", "C.Goblet", "SI.Stem", "SI.Stem", "C.Stem", "C.Goblet", 
               "C.Stem", "SI.Enterocyte", "C.Goblet", "Endocrine", "SI.Enterocyte", 
               "SI.Enterocyte", "C.Goblet", "C.Goblet", "SI.Stem", "C.Goblet", 
               "C.Stem", "C.Goblet", "SI.Stem", "SI.Enterocyte", "C.Stem", "C.Goblet", 
               "C.Goblet", "C.Stem", "SI.Stem", "C.Goblet", "C.Goblet", "Blood", 
               "C.Stem", "SI.Stem", "C.Stem", "SI.Enterocyte", "C.Stem", "C.Stem", 
               "C.Goblet", "C.Stem", "C.Goblet", "Endocrine", "SI.Tufft", "C.Colonocyte", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Goblet", "C.Stem", "C.Colonocyte", 
               "C.Colonocyte", "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", 
               "C.Goblet", "C.Stem", "C.Stem", "C.Goblet", "SI.Stem", "C.Stem", 
               "SI.Stem", "C.Colonocyte", "C.Stem", "SI.Stem", "SI.Stem", "Blood", 
               "C.Stem", "C.Colonocyte", "C.Goblet", "C.Stem", "SI.Enterocyte", 
               "C.Goblet", "SI.Stem", "C.Goblet", "C.Stem", "SI.Stem", "Endocrine", 
               "C.Goblet", "SI.Stem", "C.Goblet", "C.Goblet", "SI.Stem", "Blood", 
               "C.Stem", "C.Stem", "Endocrine", "C.Stem", "C.Stem", "SI.Stem", 
               "C.Goblet", "C.Stem", "SI.Stem", "C.Colonocyte", "SI.Goblet", 
               "C.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Goblet", "SI.Stem", 
               "C.Stem", "SI.Enterocyte", "Blood", "C.Stem", "C.Stem", "C.Goblet", 
               "SI.Enterocyte", "SI.Enterocyte", "C.Stem", "C.Goblet", "C.Stem", 
               "C.Stem", "SI.Stem", "C.Colonocyte", "C.Goblet", "C.Goblet", 
               "Endocrine", "SI.Goblet", "C.Stem", "C.Stem", "C.Goblet", "C.Goblet", 
               "Blood", "C.Stem", "C.Goblet", "Blood", "Blood", "SI.Stem", "Blood", 
               "C.Goblet", "C.Colonocyte", "Blood", "C.Stem", "C.Stem", "C.Goblet", 
               "SI.Enterocyte", "SI.Stem", "C.Colonocyte", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Enterocyte", "C.Stem", "SI.Tufft", "C.Stem", 
               "C.Stem", "C.Stem", "SI.Enterocyte", "C.Goblet", "C.Stem", "SI.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "SI.Stem", "C.Goblet", "C.Stem", 
               "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", "C.Goblet", 
               "C.Stem", "SI.Enterocyte", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "Endocrine", "C.Goblet", "C.Goblet", "C.Stem", "C.Goblet", "C.Goblet", 
               "SI.Stem", "C.Stem", "C.Stem", "C.Goblet", "C.Stem", "SI.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "C.Goblet", "Blood", "C.Goblet", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", "SI.Goblet", "C.Goblet", 
               "SI.Enterocyte", "C.Colonocyte", "C.Stem", "SI.Stem", "SI.Enterocyte", 
               "C.Goblet", "C.Stem", "C.Goblet", "SI.Stem", "SI.Stem", "C.Goblet", 
               "C.Stem", "SI.Stem", "SI.Enterocyte", "C.Stem", "SI.Stem", "SI.Enterocyte", 
               "C.Goblet", "C.Stem", "SI.Stem", "C.Colonocyte", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Enterocyte", "C.Colonocyte", "SI.Enterocyte", "SI.Enterocyte", 
               "SI.Tufft", "C.Stem", "C.Stem", "C.Colonocyte", "C.Colonocyte", 
               "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "C.Stem", "C.Stem", 
               "SI.Enterocyte", "C.Stem", "C.Goblet", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "C.Goblet", "SI.Stem", "Endocrine", "C.Colonocyte", 
               "C.Goblet", "C.Colonocyte", "SI.Stem", "C.Goblet", "C.Stem", 
               "SI.Enterocyte", "C.Goblet", "SI.Enterocyte", "SI.Stem", "C.Goblet", 
               "C.Colonocyte", "SI.Enterocyte", "C.Stem", "C.Stem", "C.Stem", 
               "SI.Stem", "SI.Enterocyte", "C.Stem", "C.Stem", "SI.Enterocyte", 
               "C.Colonocyte", "Blood", "Blood", "C.Stem", "C.Stem", "C.Stem", 
               "SI.Stem", "SI.Stem", "C.Goblet", "C.Stem", "C.Goblet", "C.Colonocyte", 
               "C.Goblet", "C.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "Endocrine", "C.Goblet", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", 
               "C.Goblet", "SI.Stem", "C.Goblet", "C.Colonocyte", "Blood", "SI.Stem", 
               "SI.Stem", "C.Stem", "C.Colonocyte", "SI.Enterocyte", "C.Stem", 
               "SI.Enterocyte", "SI.Goblet", "SI.Stem", "C.Goblet", "C.Stem", 
               "SI.Enterocyte", "C.Colonocyte", "C.Goblet", "C.Goblet", "C.Stem", 
               "C.Goblet", "Blood", "SI.Enterocyte", "C.Goblet", "Blood", "SI.Enterocyte", 
               "C.Stem", "SI.Stem", "C.Colonocyte", "C.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "Blood", "SI.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "SI.Goblet", "SI.Enterocyte", 
               "C.Colonocyte", "C.Colonocyte", "C.Stem", "Endocrine", "SI.Enterocyte", 
               "SI.Enterocyte", "C.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "C.Stem", "C.Goblet", "C.Goblet", "C.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "C.Goblet", "C.Colonocyte", "SI.Enterocyte", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Goblet", "C.Goblet", "C.Goblet", 
               "SI.Stem", "C.Colonocyte", "C.Colonocyte", "Endocrine", "C.Stem", 
               "SI.Tufft", "Endocrine", "C.Stem", "C.Stem", "Endocrine", "C.Stem", 
               "C.Goblet", "C.Stem", "C.Stem", "SI.Stem", "C.Goblet", "SI.Stem", 
               "C.Goblet", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", 
               "C.Stem", "C.Goblet", "SI.Stem", "SI.Stem", "Blood", "SI.Stem", 
               "C.Stem", "C.Goblet", "C.Goblet", "C.Goblet", "C.Colonocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "C.Stem", "C.Goblet", 
               "SI.Stem", "SI.Stem", "C.Goblet", "C.Stem", "C.Stem", "C.Stem", 
               "C.Stem", "C.Stem", "C.Colonocyte", "C.Stem", "SI.Stem", "C.Goblet", 
               "SI.Enterocyte", "C.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Goblet", 
               "SI.Enterocyte", "C.Stem", "Endocrine", "SI.Enterocyte", "C.Colonocyte", 
               "SI.Stem", "C.Colonocyte", "SI.Goblet", "C.Stem", "SI.Enterocyte", 
               "SI.Stem", "C.Goblet", "C.Goblet", "SI.Enterocyte", "C.Stem", 
               "SI.Enterocyte", "Blood", "SI.Goblet", "C.Stem", "C.Stem", "C.Colonocyte", 
               "C.Colonocyte", "Blood", "C.Stem", "C.Goblet", "C.Colonocyte", 
               "SI.Stem", "SI.Tufft", "C.Stem", "SI.Enterocyte", "C.Colonocyte", 
               "SI.Enterocyte", "SI.Tufft", "SI.Stem", "C.Stem", "C.Stem", "C.Colonocyte", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Colonocyte", "SI.Stem", 
               "C.Stem", "C.Goblet", "C.Colonocyte", "C.Stem", "C.Colonocyte", 
               "C.Stem", "C.Stem", "C.Stem", "C.Goblet", "C.Stem", "C.Goblet", 
               "C.Colonocyte", "C.Stem", "C.Colonocyte", "C.Goblet", "C.Stem", 
               "C.Stem", "SI.Stem", "Blood", "C.Stem", "C.Stem", "C.Goblet", 
               "C.Stem", "C.Stem", "C.Goblet", "SI.Stem", "C.Colonocyte", "C.Stem", 
               "SI.Enterocyte", "C.Stem", "C.Stem", "SI.Goblet", "SI.Stem", 
               "C.Stem", "Blood", "SI.Enterocyte", "C.Stem", "C.Goblet", "C.Goblet", 
               "C.Goblet", "C.Stem", "C.Colonocyte", "SI.Enterocyte", "SI.Stem", 
               "C.Colonocyte", "C.Stem", "SI.Enterocyte", "C.Stem", "C.Stem", 
               "SI.Stem", "C.Colonocyte", "C.Stem", "SI.Stem", "C.Stem", "C.Stem", 
               "C.Stem", "C.Goblet", "Blood", "SI.Enterocyte", "C.Stem", "C.Stem", 
               "SI.Enterocyte", "C.Stem", "SI.Stem", "C.Colonocyte", "SI.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "Blood", "SI.Stem", "C.Stem", 
               "C.Goblet", "SI.Stem", "SI.Tufft", "C.Stem", "C.Goblet", "C.Stem", 
               "C.Colonocyte", "C.Stem", "SI.Stem", "C.Goblet", "C.Stem", "C.Stem", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Goblet", "C.Stem", 
               "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "SI.Tufft", "C.Stem", "SI.Enterocyte", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "C.Stem", "C.Goblet", "SI.Tufft", "SI.Stem", 
               "C.Stem", "Endocrine", "C.Goblet", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "C.Stem", "SI.Stem", "C.Goblet", "C.Stem", "C.Stem", 
               "C.Stem", "SI.Stem", "C.Goblet", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Goblet", "SI.Stem", "SI.Stem", "C.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Goblet", 
               "SI.Enterocyte", "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Tufft", "C.Stem", "SI.Stem", "C.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "C.Goblet", "SI.Stem", "C.Goblet", "C.Stem", "SI.Goblet", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "Endocrine", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "Endocrine", 
               "C.Stem", "SI.Stem", "SI.Tufft", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "C.Stem", "C.Goblet", "SI.Stem", "SI.Enterocyte", 
               "SI.Goblet", "C.Stem", "SI.Stem", "C.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Enterocyte", "C.Goblet", "C.Stem", "SI.Stem", 
               "C.Goblet", "SI.Stem", "C.Goblet", "C.Stem", "C.Stem", "SI.Goblet", 
               "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "C.Goblet", 
               "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Enterocyte", "C.Stem", 
               "SI.Goblet", "C.Goblet", "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "C.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "C.Stem", "C.Goblet", 
               "C.Goblet", "C.Stem", "C.Stem", "C.Stem", "SI.Goblet", "SI.Enterocyte", 
               "SI.Stem", "C.Stem", "SI.Goblet", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "C.Stem", "SI.Stem", "SI.Tufft", "C.Stem", "C.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "C.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", "SI.Goblet", 
               "C.Stem", "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "C.Goblet", 
               "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", "C.Stem", 
               "SI.Stem", "C.Stem", "SI.Stem", "SI.Enterocyte", "C.Stem", "C.Goblet", 
               "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", 
               "C.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", "Endocrine", 
               "SI.Tufft", "C.Stem", "C.Stem", "C.Stem", "SI.Enterocyte", "C.Goblet", 
               "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", "C.Goblet", "C.Stem", 
               "SI.Tufft", "SI.Enterocyte", "C.Stem", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", 
               "SI.Enterocyte", "C.Stem", "C.Goblet", "C.Stem", "C.Stem", "C.Goblet", 
               "C.Stem", "SI.Goblet", "C.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "C.Goblet", "SI.Enterocyte", "C.Stem", "SI.Tufft", 
               "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "C.Stem", "SI.Stem", "C.Stem", "C.Stem", 
               "C.Goblet", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", "C.Stem", 
               "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Goblet", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Enterocyte", "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "C.Goblet", "C.Goblet", "SI.Enterocyte", "SI.Stem", 
               "C.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "C.Stem", "SI.Stem", "C.Stem", "SI.Enterocyte", "C.Stem", 
               "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", 
               "SI.Stem", "C.Stem", "SI.Goblet", "SI.Stem", "SI.Goblet", "C.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "Endocrine", "C.Stem", "C.Stem", 
               "C.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "C.Goblet", 
               "C.Stem", "C.Stem", "SI.Goblet", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", 
               "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Enterocyte", "C.Stem", "SI.Stem", "C.Goblet", 
               "C.Stem", "C.Stem", "C.Stem", "SI.Stem", "C.Stem", "C.Stem", 
               "SI.Stem", "Endocrine", "SI.Enterocyte", "C.Stem", "SI.Tufft", 
               "C.Stem", "SI.Stem", "C.Stem", "SI.Tufft", "C.Goblet", "C.Stem", 
               "C.Goblet", "C.Stem", "SI.Stem", "C.Goblet", "SI.Enterocyte", 
               "C.Stem", "C.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "C.Stem", 
               "C.Stem", "SI.Tufft", "SI.Stem", "SI.Enterocyte", "C.Stem", "C.Stem", 
               "C.Stem", "Endocrine", "C.Stem", "C.Goblet", "C.Goblet", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "C.Stem", 
               "SI.Goblet", "C.Goblet", "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", 
               "SI.Stem", "C.Goblet", "C.Stem", "C.Goblet", "C.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "C.Goblet", "C.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "C.Stem", 
               "C.Goblet", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", "C.Goblet", 
               "C.Stem", "Endocrine", "SI.Enterocyte", "C.Stem", "C.Stem", "SI.Goblet", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", 
               "SI.Goblet", "C.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "SI.Tufft", "C.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "C.Stem", "C.Stem", "SI.Goblet", "C.Stem", "C.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Enterocyte", "SI.Goblet", "C.Stem", 
               "C.Goblet", "C.Goblet", "SI.Goblet", "SI.Goblet", "SI.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "SI.Enterocyte", "C.Goblet", "C.Goblet", 
               "C.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", 
               "C.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", 
               "C.Stem", "C.Goblet", "C.Stem", "SI.Enterocyte", "C.Stem", "C.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "C.Stem", "SI.Enterocyte", "C.Stem", "SI.Enterocyte", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "C.Goblet", "C.Goblet", 
               "C.Stem", "SI.Stem", "SI.Stem", "C.Goblet", "C.Goblet", "C.Stem", 
               "C.Goblet", "C.Goblet", "SI.Enterocyte", "SI.Stem", "C.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Tufft", "C.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "C.Goblet", "SI.Tufft", "C.Stem", "C.Goblet", "SI.Stem", "C.Stem", 
               "SI.Goblet", "SI.Goblet", "SI.Stem", "Endocrine", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "Blood", "SI.Tufft", "SI.Stem", "C.Goblet", 
               "SI.Stem", "SI.Goblet", "C.Stem", "C.Stem", "Endocrine", "C.Stem", 
               "Endocrine", "C.Stem", "SI.Stem", "C.Stem", "C.Goblet", "SI.Stem", 
               "SI.Enterocyte", "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", 
               "Endocrine", "C.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "C.Stem", "SI.Goblet", "C.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "C.Goblet", "C.Goblet", "C.Goblet", "SI.Goblet", "C.Goblet", 
               "C.Stem", "C.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Tufft", 
               "SI.Stem", "C.Stem", "C.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "C.Goblet", "C.Stem", "C.Goblet", "C.Stem", "C.Stem", "C.Stem", 
               "SI.Goblet", "C.Stem", "SI.Stem", "C.Goblet", "SI.Stem", "SI.Stem", 
               "C.Stem", "SI.Goblet", "C.Stem", "C.Goblet", "SI.Stem", "C.Stem", 
               "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "SI.Stem", "C.Goblet", "Endocrine", "C.Goblet", "C.Goblet", "SI.Stem", 
               "SI.Goblet", "C.Stem", "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", 
               "C.Stem", "SI.Tufft", "C.Goblet", "SI.Goblet", "C.Stem", "SI.Goblet", 
               "C.Goblet", "C.Stem", "SI.Stem", "SI.Stem", "C.Stem", "SI.Goblet", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "C.Stem", "C.Goblet", "SI.Stem", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "SI.Goblet", "C.Goblet", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "C.Goblet", "C.Goblet", "SI.Stem", 
               "C.Stem", "C.Goblet", "SI.Stem", "SI.Stem", "C.Goblet", "C.Goblet", 
               "Blood", "SI.Stem", "C.Goblet", "C.Goblet", "SI.Stem", "SI.Enterocyte", 
               "Blood", "C.Goblet", "C.Goblet", "C.Goblet", "C.Goblet", "Endocrine", 
               "C.Goblet", "SI.Stem", "SI.Goblet", "C.Goblet", "C.Goblet", "C.Stem", 
               "C.Stem", "C.Stem", "SI.Stem", "C.Goblet", "C.Stem", "C.Goblet", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "C.Stem", "C.Goblet", 
               "C.Stem", "SI.Goblet", "C.Stem", "SI.Stem", "Blood", "SI.Stem", 
               "SI.Paneth", "C.Stem", "SI.Stem", "C.Goblet", "SI.Stem", "SI.Goblet", 
               "C.Stem", "SI.Goblet", "C.Stem", "C.Goblet", "C.Goblet", "SI.Stem", 
               "Blood", "C.Goblet", "C.Stem", "C.Goblet", "C.Goblet", "SI.Stem", 
               "C.Goblet", "SI.Stem", "C.Stem", "SI.Tufft", "C.Goblet", "C.Stem", 
               "C.Stem", "SI.Stem", "C.Stem", "C.Stem", "C.Stem", "SI.Stem", 
               "C.Stem", "SI.Goblet", "SI.Tufft", "SI.Goblet", "C.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Enterocyte", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Stem", "SI.Enterocyte", "SI.Paneth", "SI.Enterocyte", "SI.Tufft", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Tufft", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Paneth", "SI.Enterocyte", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", 
               "SI.Goblet", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Enterocyte", "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Paneth", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Tufft", 
               "SI.Paneth", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Goblet", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "Endocrine", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Paneth", "SI.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Paneth", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Stem", 
               "Endocrine", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Tufft", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "Endocrine", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Enterocyte", 
               "SI.Tufft", "SI.Goblet", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Enterocyte", "Endocrine", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Paneth", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Tufft", "SI.Paneth", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "Endocrine", "SI.Goblet", "SI.Stem", "SI.Goblet", 
               "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Paneth", "SI.Stem", "SI.Enterocyte", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", "SI.Tufft", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Paneth", "SI.Enterocyte", "SI.Paneth", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Tufft", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "C.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Paneth", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "Endocrine", "SI.Tufft", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Paneth", "SI.Stem", "SI.Stem", "Endocrine", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Goblet", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Goblet", "SI.Enterocyte", 
               "SI.Goblet", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "Endocrine", 
               "SI.Stem", "SI.Goblet", "SI.Enterocyte", "SI.Stem", "SI.Paneth", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Paneth", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Enterocyte", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Enterocyte", 
               "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "Endocrine", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Paneth", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", "SI.Stem", 
               "SI.Paneth", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Enterocyte", 
               "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Paneth", "SI.Tufft", "SI.Stem", "C.Stem", "SI.Enterocyte", 
               "SI.Enterocyte", "SI.Stem", "SI.Goblet", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Goblet", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "C.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "Endocrine", "SI.Goblet", "SI.Stem", "SI.Goblet", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", "SI.Stem", "SI.Paneth", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Paneth", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", 
               "SI.Paneth", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Goblet", "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Goblet", 
               "SI.Stem", "SI.Paneth", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Tufft", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "Endocrine", "SI.Paneth", "SI.Enterocyte", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "Endocrine", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Paneth", "SI.Enterocyte", "SI.Stem", "SI.Goblet", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "Endocrine", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Paneth", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Goblet", "SI.Goblet", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Enterocyte", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Tufft", 
               "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "Endocrine", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Tufft", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Paneth", "SI.Enterocyte", 
               "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Stem", 
               "SI.Enterocyte", "SI.Paneth", "SI.Stem", "Endocrine", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Goblet", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Paneth", "SI.Paneth", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", "SI.Tufft", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Paneth", "SI.Goblet", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Paneth", "SI.Tufft", "SI.Paneth", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "Endocrine", "SI.Stem", 
               "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Stem", "SI.Paneth", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Stem", "SI.Enterocyte", "SI.Enterocyte", "SI.Stem", "SI.Enterocyte", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Enterocyte", "SI.Tufft", "SI.Stem", "SI.Tufft", 
               "SI.Enterocyte", "SI.Stem", "SI.Paneth", "SI.Goblet", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "SI.Stem", "SI.Tufft", "SI.Tufft", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "C.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Goblet", "SI.Goblet", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Stem", "SI.Enterocyte", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", "SI.Stem", 
               "Endocrine", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Goblet", "SI.Enterocyte", "SI.Stem", "SI.Stem", "SI.Stem", 
               "SI.Stem", "SI.Tufft", "SI.Paneth", "SI.Enterocyte", "SI.Stem", 
               "SI.Stem", "SI.Stem", "Endocrine", "SI.Stem", "SI.Stem", "SI.Tufft", 
               "SI.Tufft")
  
  uClasses <- sort(unique(classes))
  i <- case_when(
    classes == uClasses[1] ~ 1,
    classes == uClasses[2] ~ 2,
    classes == uClasses[3] ~ 3,
    classes == uClasses[4] ~ 4,
    classes == uClasses[5] ~ 5,
    classes == uClasses[6] ~ 6,
    classes == uClasses[7] ~ 7,
    classes == uClasses[8] ~ 8,
    classes == uClasses[9] ~ 9,
    classes == uClasses[10] ~ 10,
    TRUE ~ 0
  )
  names(i) <- uClasses
  singlets <- rbind(i, i, i, i, i)
  
  #
  output <- sp.scRNAseq:::.subsetSinglets(classes, singlets, 10)
})

c <- tibble(
  sample = rownames(getData(uObj, "tsne")),
  class = getData(uObj, "classification")
)

getData(cObjSng, "counts.cpm") %>%
  matrix_to_tibble("gene") %>%
  filter(gene == "Lgr5") %>%
  gather(sample, value, -gene) %>%
  inner_join(c) %>%
  ggplot() +
  geom_boxplot(aes(class, value))

sp.scRNAseq:::.subsetSinglets(classes, singlets, 400) %>%
  matrix_to_tibble("gene") %>%
  filter(gene == "Lgr5") %>%
  gather(class, value, -gene) %>%
  ggplot() +
  geom_point(aes(class, value))

