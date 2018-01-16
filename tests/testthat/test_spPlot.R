#context("spPlot")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[ ,s], testErcc[ ,s])
cObjMul <- spCounts(testCounts[ ,!s], testErcc[ ,!s])
uObj <- testUns
sObj <- testSwa

################################################################################
#                                                                              #
# Counts Plots                                                                 #
#                                                                              #
################################################################################

##run test .countsErccPlot
test_that("check that .countsErccPlot outputs the expected result", {
    
    ###TEST1####
    #run function
    output <- .countsErccPlot(cObjSng, cObjMul)
    
    #test
    expect_is(output$layers[[1]], "ggproto")
    expect_identical(
        class(output$layers[[1]]$geom),
        c("GeomPoint","Geom","ggproto")
    )
    expect_identical(output$labels$x, "Sample type")
    expect_identical(output$labels$y, "Cell number")
    expect_identical(output$labels$title, "Fraction ERCC in singlets/doublets")
    expect_false(is.null(output))
    expect_silent(.countsErccPlot(cObjSng, cObjMul))

    ###TEST2###
    #here you should probably try to trigger this monotome 2nd axis error from
    #ggplot to see when it arises
    
})

##run test .countsMarkersPlotProcess
test_that("check that .countsMarkersPlotProcess outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    markers <- c("a1", "b1")
    
    #setup expected data
    expected <- data.frame(
        sampleType = c(rep("Singlet", 340), rep("Multuplet", 2)),
        marker1 = c(
            getData(cObjSng, "counts.log")['a1',],
            getData(cObjMul, "counts.log")['a1',]
        ),
        marker2 = c(
            getData(cObjSng, "counts.log")['b1',],
            getData(cObjMul, "counts.log")['b1',]
        )
    )
    
    #run function
    output <- .countsMarkersPlotProcess(cObjSng, cObjMul, markers)
    
    #test
    expect_equivalent(expected, output)
    expect_type(output, "list")
    expect_type(output$sampleType, "integer")
    expect_type(output$marker1, "double")
    expect_type(output$marker2, "double")
})

##run test .countsMarkersPlot
test_that("check that .countsMarkersPlot outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    markers <- c("a1", "b1")
    
    #run function
    output <- .countsMarkersPlot(cObjSng, cObjMul, markers)
    
    #test
    expect_is(output$layers[[1]], "ggproto")
    expect_identical(
        class(output$layers[[1]]$geom),
        c("GeomPoint","Geom","ggproto")
    )
    expect_identical(output$labels$x, "log2( Normalized counts: a1 )")
    expect_identical(output$labels$y, "log2( Normalized counts: b1 )")
    expect_identical(output$labels$title, "Cell Identity Markers")
    expect_identical(output$labels$colour, "sampleType")
    expect_identical(output$scales$scales[[1]]$name, "sampleType")
    expect_false(is.null(output))
    expect_silent(.countsMarkersPlot(cObjSng, cObjMul, markers))
})

################################################################################
#                                                                              #
# Unsupervised Plots                                                           #
#                                                                              #
################################################################################

##run test .unsupClusterPlotProcess
test_that("check that .unsupClusterPlotProcess outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- uObj
    
    #run function
    output <- .unsupClusterPlotProcess(input, plotUncertainty=FALSE)
    
    #test
    expect_equivalent(nrow(output), 340)
    expect_equivalent(ncol(output), 5)
    expect_identical(
        colnames(output),
        c("sample", "V1", "V2", "classification", "uncertainty")
    )
    expect_type(output, "list")
    expect_type(output$V1, "double")
    expect_type(output$V2, "double")
    expect_type(output$classification, "character")
    expect_type(output$uncertainty, "double")
})

##run test .unsupClustersPlot
test_that("check that .unsupClustersPlot outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- uObj

    #run function
    output <- .unsupClustersPlot(input, plotUncertainty=FALSE)
    
    #test
    expect_identical(output$labels$x, "Dim 1")
    expect_identical(output$labels$y, "Dim 2")
    expect_identical(output$labels$title, "Clusters")
    expect_identical(output$labels$colour, "classification")
    expect_identical(output$scales$scales[[1]]$name, "sampleType")
    expect_silent(.unsupClustersPlot(input, plotUncertainty=FALSE))
    expect_false(is.null(output))
    
})

##run test .unsupMarkerPlotProcess
test_that("check that .unsupMarkerPlotProcess outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    markers <- c("a1", "b1")
    
    #run function
    output <- .unsupMarkerPlotProcess(
        uObj,
        cObjSng,
        markers,
        plotUncertainty=FALSE
    )
    
    #test
    expect_equivalent(nrow(output), 340*2)
    expect_equivalent(ncol(output), 6)
    cols <- c("sample", "uncertainty", "V1", "V2", "variable", "value")
    expect_identical(colnames(output), cols)
    expect_type(output$V1, "double")
    expect_type(output$V2, "double")
    expect_type(output$sample, "character")
    expect_type(output$uncertainty, "double")
    expect_type(output$variable, "character")
    expect_type(output$value, "double")

})

##run test .unsupMarkersPlot
test_that("check that the .unsupMarkersPlot function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    markers <- c("a1", "b1")

    #run function
    output <- .unsupMarkersPlot(
        uObj,
        cObjSng,
        markers,
        plotUncertainty=FALSE
    )
    
    #test
    expect_identical(output$labels$x, "x")
    expect_identical(output$labels$y, "y")
    expect_identical(output$labels$title, "Markers")
    expect_identical(output$labels$colour, "Expression")
    expect_silent(.unsupMarkersPlot(uObj,cObjSng,markers,plotUncertainty=FALSE))
    expect_false(is.null(output))
    
})

################################################################################
#                                                                              #
# Swarm Plots                                                                  #
#                                                                              #
################################################################################
