#context("spSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])
uObj <- testUns
sObj <- testSwa

##run test .defineImport
test_that("check that    outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    
    #setup expected data
    
    #run function

    #test

})

##run test getMultipletsForEdge
test_that("check that getMultipletsForEdge outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    #A1 and B1 should have an edge
    #I1 and J1 should have an edge
    expected1 <- c("A1-B1" = "m.A1B1")
    expected2 <- c("C1-D1" = "m.C1D1")
    expected3 <- c("A1-B1" = "m.A1B1", "C1-D1" = "m.C1D1")
    
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
    expected <- list(
        "A1-B1" = c("m.A1B1C1a", "m.A1B1C1b"),
        "A1-C1" = "m.A1B1C1a",
        "A1-D1" = "m.A1B1C1b",
        "B1-C1" = "m.A1B1C1a",
        "B1-D1" = "m.A1B1C1b",
        "C1-D1" = character()
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
    expected1 <- data.frame(from="A1", to="B1", stringsAsFactors=FALSE)
    expected2 <- data.frame(from="C1", to="D1", stringsAsFactors=FALSE)

    #run function
    output1 <- getEdgesForMultiplet(sObj, 1/4, 'm.A1B1')[,1:2]
    output2 <- getEdgesForMultiplet(sObj, 1/4, 'm.C1D1')[,1:2]
    rownames(output1) <- NULL
    rownames(output2) <- NULL

    #test
    expect_identical(output1, expected1)
    expect_identical(output2, expected2)

})

##run test .makePermutations
test_that("check that .makePermutations outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    classes <- LETTERS[1:5]
    iter <- 10
    
    #setup expected data
    
    #run function
    output <- .makePermutations(classes, iter, seed=11)
    
    #test
    expect_true(nrow(output) == length(LETTERS[1:5]))
    expect_true(ncol(output) == iter)
    expect_true(class(output) == "matrix")
})

##run test .runPermutations
test_that("check that .runPermutations outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    iter <- 5
    classes <- getData(testUns, "classification")
    maxiter <- 5
    swarmsize <- 100
    cores <- 1
    seed <- 11
    norm <- TRUE
    permMatrix <- .makePermutations(classes, iter, seed)
    
    #setup expected data
    
    #run function
    output <- .runPermutations(
        permMatrix = permMatrix,
        iter = iter,
        spCounts = cObjSng,
        spCountsMul = cObjMul,
        spUnsupervised = testUns,
        classes = classes,
        distFun = distToSlice,
        maxiter = maxiter,
        swarmsize = swarmsize,
        cores = cores,
        seed = seed,
        norm = norm
    )
    
    #test
    expect_identical(dim(output), c(10L,5L))
    expect_true(class(output) == "matrix")
    expect_type(output, "integer")
})

##run test .calculatePermP
test_that("check that .calculatePermP outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    iter <- 10^4
    permData <- matrix(
        c(
            c(rep(0, 10^4)), #NA
            c(rep(0, 5*10^3), rep(1, 5*10^3)), #NA
            c(rep(100, 10^4)), #NA
            c(rep(2, 10^4)), #NA
            c(2, rep(0, 9999)), #1/10^4
            c(1, rep(0, 9999)), #1/10^4
            c(rep(0, 10^4)), #10^-log10(iter)
            c(rep(0, 10^4)), #10^-log10(iter)
            c(rep(2, 10^4)), #1
            c(rep(0, 5*10^3), rep(1, 5*10^3)) #0.5
        ),
        byrow=TRUE,
        ncol=10^4
    )
    
    #setup expected data
    expected <- c(rep(NA, 4), rep(1/10^4, 2), rep(10^-log10(iter), 2), 1, 0.5)
    
    #run function
    output <- .calculatePermP(permData, sObj, iter)
    
    #test
    expect_identical(output$pValue, expected)
    expect_true(class(output) == "data.frame")
    expect_type(output$from, "character")
    expect_type(output$to, "character")
    expect_type(output$weight, "integer")
    expect_type(output$pValue, "double")

})


