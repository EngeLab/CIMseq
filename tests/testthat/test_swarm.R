context("spSwarm")

#Function to check if all elements in a vector are identical
has_zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

#run test getMultipletsForEdge
test_that("check that getMultipletsForEdge outputs the expected result", {

  ###TEST1####
  #setup expected data
  #A1 and B1 should have an edge
  #I1 and J1 should have an edge
  expected1 <- tibble::tibble(
    sample = "m.NJB00204.G04",
    from = "A375",
    to = "HOS"
  )
  expected2 <- tibble::tibble(
    sample = "m.NJB00204.D07",
    from = "HCT116",
    to = "HOS"
  )
  expected3 <- tibble::tibble(
    sample = c("m.NJB00204.A02", "m.NJB00204.G04"),
    from = c("A375", "A375"),
    to = c("HCT116", "HOS")
  )

    #run function
  output1 <- getMultipletsForEdge(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 
    data.frame("A375", "HOS")
  )
  output2 <- getMultipletsForEdge(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 
    data.frame("HCT116", "HOS")
  )
  output3 <- getMultipletsForEdge(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test,
    data.frame(c("A375", "A375"), c("HCT116", "HOS"))
  )

    #test
  expect_identical(output1, expected1)
  expect_identical(output2, expected2)
  expect_identical(output3, expected3)
})

##run test getEdgesForMultiplet
test_that("check that getEdgesForMultiplet outputs the expected result", {

  ###TEST1####
  #setup expected data
  expected1 <- tibble::tibble(
    sample = c("m.NJB00204.D07", "m.NJB00204.D07"),
    from =  c("HOS", "HCT116"),
    to = c("HCT116", "HOS")
  )
  expected2 <- tibble::tibble(
    sample = c("m.NJB00204.G04", "m.NJB00204.G04"),
    from = c("HOS", "A375"),
    to = c("A375", "HOS")
  )
  
  #run function
  output1 <- getEdgesForMultiplet(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 'm.NJB00204.D07'
  )
  output2 <- getEdgesForMultiplet(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 'm.NJB00204.G04'
  )

  #test
  expect_identical(output1, expected1)
  expect_identical(output2, expected2)
  
  ###TEST2####
  #setup expected data
  expected1 <- tibble::tibble(
    sample = rep(c("m.NJB00204.G04", "m.NJB00204.D07"), each = 2),
    from = c("HOS", "A375", "HOS", "HCT116"), 
    to = c("A375", "HOS", "HCT116", "HOS")
  )
  
  #run function
  output1 <- getEdgesForMultiplet(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 
    c('m.NJB00204.G04', 'm.NJB00204.D07')
  )
  
  #test
  expect_identical(output1, expected1)
})

test_that("check that adjustFractions outputs the expected result", {
  
  ###TEST1####
  #setup expected data
  fractions <- getData(CIMseqSwarm_test, "fractions")
  
  cnc <- cellNumberPerClass(CIMseqSinglets_test, CIMseqMultiplets_test) %>%
  {setNames(pull(., medianCellNumber), pull(., class))}
  cnc <- cnc[match(colnames(fractions), names(cnc))]
  
  cnm <- estimateCells(CIMseqSinglets_test, CIMseqMultiplets_test) %>%
    filter(sampleType == "Multiplet") %>%
    {setNames(pull(., estimatedCellNumber), pull(., sample))}
  cnm <- cnm[match(names(cnm), rownames(fractions))]
  
  frac.renorm <- t(t(fractions) / cnc)
  expected1 <- round(frac.renorm * cnm)
  expected2 <- expected1
  expected1[expected1 > 0] <- 1
  
  #run function
  output1 <- adjustFractions(
    CIMseqSinglets_test, CIMseqMultiplets_test, CIMseqSwarm_test, binary = TRUE
  )
  
  output2 <- adjustFractions(
    CIMseqSinglets_test, CIMseqMultiplets_test, CIMseqSwarm_test, binary = FALSE
  )
  
  #test
  expect_identical(expected1, output1)
  expect_identical(expected2, output2)
})

test_that("check that calculateEdgeStats outputs the expected result", {
  
  ###TEST1####
  #setup input data
  set.seed(93223)
  mat <- matrix(
    sample(c(0, 1), 30, replace = TRUE), ncol = 3, 
    dimnames = list(NULL, LETTERS[1:3])
  )
  
  #setup expected data
  expected <- c(2L, 5L, 2L, 3L, 5L, 3L)
  
  #run function
  output <- .calculateWeight(mat)
  
  #test
  expect_identical(expected, output$weight)
  
  ###TEST2####
  #setup input
  set.seed(9322)
  edges <- expand.grid(
    from = LETTERS[1:3], to = LETTERS[1:3], stringsAsFactors = FALSE
  )
  edges <- edges[edges[, 1] != edges[, 2], ]
  edges$weight <- sample(1:5, nrow(edges), replace = TRUE)
  mat <- matrix(
    sample(c(0, 1), 30, replace = TRUE), ncol = 3, 
    dimnames = list(NULL, LETTERS[1:3])
  )
  
  #setup expected data
  expected <- c(7, 1, 8, 2, 1, 1)
  
  #run function
  output <- CIMseq:::.calculateP(edges, mat)
  
  #test
  expect_identical(expected, round(output$expected.edges))
  expect_identical(output$score, output$weight / output$expected.edges)
})

################################################################################
#                                                                              #
#                                C++ functions                                 #
#                                                                              #
################################################################################

context("calculateCost")

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
  expect_true(!identical(sampleSinglets(classes), sampleSinglets(classes)))
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

##run test calculateCostDensity
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

##run test calculateLogRowMeans
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

##run test fixNegInf
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

##run test appropriateSinglets
test_that("check that appropriateSinglets outputs the expected result", {
  
  ###TEST1####
  #extract needed variables
  selectInd <- getData(CIMseqSwarm_test, "arguments")$features[[1]]
  singlets <- getData(CIMseqSinglets_test, "counts.cpm")
  n <- getData(CIMseqSwarm_test, "arguments")$nSyntheticMultiplets
  idx <- getData(CIMseqSwarm_test, "singletIdx")
  
  #setup expected data
  nc <- length(unique(getData(CIMseqSinglets_test, "classification")))
  nr <- prod(length(selectInd), n)
  rn <- paste(rep(rownames(singlets)[selectInd], each = n), 1:n, sep = ".")
  singlets <- singlets[selectInd, ]
  
  col.first.first.gene <- singlets[1, idx[[1]] + 1]
  col.last.first.gene <- singlets[1, idx[[length(idx)]] + 1]
  col.first.last.gene <- singlets[nrow(singlets), idx[[1]] + 1]
  col.last.last.gene <- singlets[nrow(singlets), idx[[length(idx)]] + 1]
  
  #run function
  set.seed(82390)
  output <- appropriateSinglets(CIMseqSinglets_test, idx, selectInd)
  
  #test
  expect_equal(nr, nrow(output))
  expect_equal(nc, ncol(output))
  expect_identical(rn, rownames(output))
  expect_identical(unname(col.first.first.gene), unname(output[1, ]))
  expect_identical(unname(col.last.first.gene), unname(output[400, ]))
  expect_identical(unname(col.first.last.gene), unname(output[nrow(output) - 399, ]))
  expect_identical(unname(col.last.last.gene), unname(output[nrow(output), ]))
  expect_true(all(colnames(output) == sort(colnames(output))))
})

test_that("check that calculateCost and cost give identical results", {
  
  ###TEST1####
  #prepare normal input data
  selectInd <- getData(CIMseqSwarm_test, "arguments")$features[[1]]
  singletIdx <- getData(CIMseqSwarm_test, "singletIdx")
  singletSubset <- appropriateSinglets(CIMseqSinglets_test, singletIdx, selectInd)
  
  m <- "m.NJB00204.G04"
  oneMultiplet <- ceiling(getData(CIMseqMultiplets_test, "counts.cpm")[, m])
  fractions <- as.numeric(getData(CIMseqSwarm_test, "fractions")[m, ])
  n <- getData(CIMseqSwarm_test, "arguments")$nSyntheticMultiplets
  
  adj <- adjustAccordingToFractions(fractions, singletSubset)
  rm <- multipletSums(adj)
  sm <- vecToMat(rm, length(oneMultiplet), n)
  
  #expected
  expected <- cost(oneMultiplet, singletSubset, fractions, n)
  expected2 <- .costCalculationR(oneMultiplet, sm)
  
  #output
  output <- calculateCost(oneMultiplet, singletSubset, fractions, n)
  
  #test
  expect_equal(expected, output)
  expect_equal(expected2, output)
  
  #check with "simple example"
  #setup input data
  classes <- rep(letters[1:2], each = 2)
  singlets <- matrix(c(rep(1:10, 2), rep(11:20, 2)), ncol = 4)
  rownames(singlets) <- LETTERS[1:nrow(singlets)]
  n <- 2
  
  sm <- purrr::map(1:n, ~sampleSinglets(classes)) %>%
    purrr::map(., ~subsetSinglets(singlets, .x)) %>%
    purrr::map2(., 1:n, function(x, i) {
      rownames(x) <- paste(rownames(singlets), i, sep = ".")
      x
      }) %>%
    purrr::map(., function(x) {colnames(x) <- sort(unique(classes)); x}) %>%
    do.call("rbind", .) %>%
    .[order(rownames(.)), ]
  
  fractions <- rep(0.5, 2)
  multiplets <- matrix(rep(rowSums(singlets[,2:3]), 2), ncol = 2)
  oneMultiplet <- multiplets[, 1]
  
  #calculate cost "by hand"
  adj <- t(t(sm) * fractions)
  rs <- rowSums(adj)
  sm.hand <- matrix(rs, ncol = 2, byrow = TRUE)
  dp <- dpois(round(oneMultiplet), sm.hand)
  rm <- rowMeans(dp)
  lrm <- log10(rm)
  lrm[is.infinite(lrm)] <- -323.0052
  cost <- sum(lrm) * (-1)
  
  #calculate cost with c++ fun
  cost.cpp <- calculateCost(oneMultiplet, sm, fractions, n)
  
  #test
  expect_equal(cost.cpp, cost)
})
