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
    
    ###TEST1####
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