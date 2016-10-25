#context("spPlot")

cObj <- spCounts(testData, matrix(), "m.")
uObj <- spUnsupervised(cObj, max=250, max_iter=1000)
sObj <- spSwarm(uObj, swarmsize = 150, cores=1, cutoff=0.14)

##run test .countsErccPlotProcess
test_that("check that the .countsErccPlotProcess function outputs the expected result", {
    
    #note this test is not great at the moment since the example data has no ercc.counts
    
    ###TEST1####
    #prepare normal input data
    input <- cObj
    
    #setup expected data
    expected <- data.frame(
        sampleType = c(rep("Singlet", 850), rep("Multuplet", 2)),
        frac.ercc = as.numeric(rep(NA, 852))
    )
    
    #run function
    output <- .countsErccPlotProcess(input)
    
    #test
    expect_identical(expected, output)

})

##run test .countsErccPlot
test_that("check that the .countsErccPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- cObj

    #run function
    output <- .countsErccPlot(input)
    
    #test
    expect_silent(.countsErccPlot(input))
    expect_false(is.null(output))
    
})

##run test .countsMarkersPlotProcess
test_that("check that the .countsMarkersPlotProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- cObj
    markers <- c("a1", "b1")
    
    #setup expected data
    expected <- data.frame(
        sampleType = c(rep("Singlet", 850), rep("Multuplet", 2)),
        marker1 = getData(cObj, "counts.log")['a1',],
        marker2 = getData(cObj, "counts.log")['b1',]
    )
    
    #run function
    output <- .countsMarkersPlotProcess(input, markers)
    
    #test
    expect_equivalent(expected, output)
    
})

##run test .countsMarkersPlot
test_that("check that the .countsMarkersPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- cObj
    markers <- c("A", "B")
    
    #run function
    output <- .countsMarkersPlot(input, markers)
    
    #test
    expect_warning(.countsMarkersPlot(input, markers), regexp = NA)
    expect_false(is.null(output))
    
})

##run test .unsupClusterPlotProcess
test_that("check that the .unsupClusterPlotProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- uObj
    
    #run function
    output <- .unsupClusterPlotProcess(input)
    
    #test
    expect_equivalent(nrow(output), 4)
    expect_equivalent(ncol(output), 3)
    expect_equivalent(colnames(output), c("V1", "V2", "classification"))
    expect_type(output$V1, "double")
    expect_type(output$V2, "double")
    expect_type(output$classification, "integer")

})

##run test .unsupClustersPlot
test_that("check that the .unsupClustersPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- uObj

    #run function
    output <- .unsupClustersPlot(input)
    
    #test
    expect_warning(.unsupClustersPlot(input), regexp = NA)
    expect_false(is.null(output))
    
})

##run test .unsupMarkerPlotProcess
test_that("check that the .unsupMarkerPlotProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- uObj
    markers <- c("a1", "b1")
    
    #run function
    output <- .unsupMarkerPlotProcess(input, markers)
    
    #test
    expect_equivalent(nrow(output), 852*2)
    expect_equivalent(ncol(output), 5)
    expect_equivalent(colnames(output), c("V1", "V2", "sample", "variable", "value"))
    expect_type(output$V1, "double")
    expect_type(output$V2, "double")
    expect_type(output$sample, "character")
    expect_type(output$variable, "integer")
    expect_type(output$value, "double")

})

##run test .unsupMarkersPlot
test_that("check that the .unsupMarkersPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- uObj
    markers <- c("A", "B")

    #run function
    output <- .unsupMarkersPlot(input, markers)
    
    #test
    expect_warning(.unsupMarkersPlot(input, markers), regexp = NA)
    expect_false(is.null(output))
    
})

##run test .swarmTsneMeansProcess
test_that("check that the .swarmTsneMeansProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- nObj
    
    #run function
    output <- .swarmTsneMeansProcess(input)
    
    #test
    expect_equivalent(nrow(output), 2)
    expect_equivalent(ncol(output), 3)
    expect_equivalent(colnames(output), c("name", "x", "y"))
    expect_type(output$name, "character")
    expect_type(output$x, "double")
    expect_type(output$y, "double")
    
})

##run test .swarmTsneProcess
test_that("check that the .swarmTsneProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- nObj
    
    #run function
    output <- .swarmTsneProcess(input)
    
    #test
    expect_warning(.swarmTsneProcess(input), regexp = NA)
    expect_false(is.null(output))
    
})

##run test .swarmNetworkDF
test_that("check that the .swarmNetworkDF function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    codedSwarm <- data.frame(
        fopt = 0,
        sampleName = "m.1",
        A1 = 0,
        B1 = 0,
        C1 = 0,
        D1 = 0.5,
        E1 = 0.5
    )

    #setup expected data
    expected <- data.frame(
        from = "D1",
        to = "E1",
        stringsAsFactors=FALSE
    )
    
    #run function
    output <- .swarmNetworkDF(codedSwarm)
    
    #test
    expect_equivalent(output, expected)
    
    ###TEST2####
    #prepare normal input data
    codedSwarm <- data.frame(
        fopt = 0,
        sampleName = "m.1",
        A1 = 1/5,
        B1 = 1/5,
        C1 = 1/5,
        D1 = 1/5,
        E1 = 1/5
    )
    
    #setup expected data
    tmp <- combn(paste(LETTERS[1:5], 1, sep=""), 2)
    expected <- data.frame(
        from = tmp[1,],
        to = tmp[2,],
        stringsAsFactors=FALSE
    )
    
    #run function
    output <- .swarmNetworkDF(codedSwarm)
    
    #test
    expect_equivalent(output, expected)

})

##run test .swarmCalculateWeights
test_that("check that the .swarmCalculateWeights function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- data.frame(
        from = c("A1", "A1", "A1", "A1"),
        to = c("B1", "B1", "C1", "D1"),
        stringsAsFactors=FALSE
    )
    
    #setup expected data
    expected <- data.frame(
        from = rep("A1", 3),
        to = c("B1", "C1", "D1"),
        weight = c(2, 1, 1),
        stringsAsFactors=FALSE
    )
    expected$weight <- as.factor(expected$weight)
    
    #run function
    output <- .swarmCalculateWeights(input)
    
    #test
    expect_equivalent(output, expected)
    
})

##run test .swarmSquareTable
test_that("check that the .swarmSquareTable function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- data.frame(
        from = c("A1", "A1", "A1", "A1"),
        to = c("B1", "B1", "C1", "D1"),
        stringsAsFactors=FALSE
    )
    
    #run function
    output <- .swarmSquareTable(input$from, input$to)
    
    #test
    expect_equivalent(
        matrix(output[lower.tri(output)], ncol=3),
        matrix(output[upper.tri(output)], ncol=3, byrow=TRUE)
    )
    
})





