#context("spSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])
uObj <- testUns
sObj <- testSwa

##run test .makeSyntheticSlice
test_that("check that .makeSyntheticSlice outputs the expected result", {
  
    ###TEST1####
    #prepare normal input data
    cellTypes <- matrix(rep(1, 40), ncol = 4)
    fractions <- c(1, 0.5, 0.25, 0.1)
    
    #setup expected data
    expected <- rep(1.85, 10)
    
    #run function
    output <- .makeSyntheticSlice(cellTypes, fractions)
    
    #test
    expect_identical(output, expected)
})

##run test .optim.fn
test_that("check that .optim.fn outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  i <- 1
  fractions <- rep(1, 4)
  distFun <- distToSlice
  cellTypes <- matrix(rep(1, 40), ncol = 4)
  control <- list(maxit = 2, s = 10)
  multiplets <- matrix(rep(1, 10), ncol = 1)
  
  #setup expected data
  expected <- rep(1, 4)
  
  #run function
  output <- .optim.fn(i, fractions, distFun, cellTypes, control, multiplets, cellNumbers = NULL, e = NULL)$par
  
  #test
  expect_identical(output, expected)
})

##run test .runPyRSwarm
test_that("check that .runPyRSwarm outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  cellTypes <- matrix(rep(1, 40), ncol = 4)
  multiplets <- matrix(rep(1, 10), ncol = 1)
  fractions <- rep(1, 4)
  distFun <- distToSlice
  maxiter <- 2
  swarmsize <- 10
  cores <- 1
  seed <- 1134
  norm <- TRUE
  report <- FALSE
  cellNumbers <- NULL
  e <- NULL
  
  #setup expected data
  dat <- data.frame(0.25, 0.25, 0.25, 0.25)
  colnames(dat) <- NULL
  
  expected <- list(
    dat,
    0,
    "Maximal number of iterations reached.",
    list()
  )
  
  #run function
  output <- .runPyRSwarm(
    cellTypes,
    multiplets,
    fractions,
    distFun,
    maxiter,
    swarmsize,
    cores,
    seed,
    norm,
    report,
    cellNumbers,
    e
  )
  
  #test
  expect_identical(output, expected)
})

##run test distToSlice
test_that("check that distToSlice outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  fractions <- rep(1, 4)
  cellTypes <- matrix(rep(1, 40), ncol = 4)
  oneMultiplet <- matrix(rep(1, 10), ncol = 1)
  
  #setup expected data
  expected <- 0
  
  #run function
  output <- distToSlice(fractions, cellTypes, oneMultiplet)
  
  #test
  expect_identical(output, expected)
})

##run test distToSliceNorm
test_that("check that distToSliceNorm outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  fractions <- rep(1, 4)
  cellTypes <- matrix(rep(1, 40), ncol = 4)
  oneMultiplet <- matrix(rep(1, 10), ncol = 1)
  
  #setup expected data
  expected <- 0
  
  #run function
  output <- distToSliceNorm(fractions, cellTypes, oneMultiplet)
  
  #test
  expect_identical(output, expected)
})

##run test getMultipletsForEdge
test_that("check that getMultipletsForEdge outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    #A1 and B1 should have an edge
    #I1 and J1 should have an edge
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
    expected3 <- tibble::tibble(
        multiplet = c("m.A1B1", "m.C1D1"),
        from = c("A1", "C1"),
        to = c("B1", "D1")
    )
    
    #run function
    output1 <- getMultipletsForEdge(sObj, 1/4, data.frame("A1", "B1"))
    output2 <- getMultipletsForEdge(sObj, 1/4, data.frame("C1", "D1"))
    output3 <- getMultipletsForEdge(
        sObj,
        1/4,
        data.frame(
            c("A1", "C1"),
            c("B1", "D1")
        )
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
