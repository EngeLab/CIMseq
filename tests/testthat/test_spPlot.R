#context("spPlot")

##run test .erccPlotProcess
test_that("check that the .erccPlotProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    counts <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    counts.ercc <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    names <- c("s.1", "m.1")
    colnames(counts) <- names
    colnames(counts.ercc) <- names

    cObj <- spCounts(counts, counts.ercc, "m.")
    
    #setup expected data
    expected <- data.frame(
        sampleType = c("Singlet", "Multuplet"),
        frac.ercc = rep(0.5, 2),
        row.names = names
    )
    
    #run function
    output <- .erccPlotProcess(cObj)
    
    #test
    expect_identical(expected, output)

})

##run test .erccPlot
test_that("check that the .erccPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    counts <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    counts.ercc <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    names <- c("s.1", "m.1")
    colnames(counts) <- names
    colnames(counts.ercc) <- names
    
    cObj <- spCounts(counts, counts.ercc, "m.")
    
    #run function
    output <- .erccPlotProcess(cObj)
    
    #test
    expect_silent(.erccPlotProcess(cObj))
    expect_false(is.null(output))
    
})

##run test .markersPlotProcess
test_that("check that the .markersPlotProcess function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    counts <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    names <- c("s.1", "m.1")
    colnames(counts) <- names
    rownames(counts) <- c("A", "B")
    cObj <- spCounts(counts, matrix(), "m.")
    
    markers <- c("A", "B")
    
    #setup expected data
    expected <- data.frame(
        sampleType = c("Singlet", "Multuplet"),
        marker1 = c(0, 19.93157),
        marker2 = c(19.93157, 0),
        row.names = names
    )
    
    #run function
    output <- .markersPlotProcess(cObj, markers)
    
    #test
    expect_equivalent(expected, output)
    
})

##run test .markersPlot
test_that("check that the .markersPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    counts <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    names <- c("s.1", "m.1")
    rownames(counts) <- c("A", "B")
    colnames(counts) <- names
    
    cObj <- spCounts(counts, matrix(), "m.")
    markers <- c("A", "B")

    #run function
    output <- .markersPlot(cObj, markers)
    
    #test
    expect_warning(.markersPlot(cObj, markers), regexp = NA)
    expect_false(is.null(output))
    
})
