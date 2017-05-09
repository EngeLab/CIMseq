#context("spUnsupervised")

##run test spTopVar
test_that("check that the spTopVar function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    s <- grepl("^s", colnames(testCounts))
    cObjSng <- spCounts(testCounts[,s], testErcc[,s])
    
    #setup expected data
    expected1 <- 27L
    expected2 <- 243L
    
    #run function
    output <- spTopVar(cObjSng, 250)
    
    #test
    expect_identical(expected1, output[1])
    expect_identical(expected2, output[length(output)])
    
})

##run test spTopMax
test_that("check that the spTopMax function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    s <- grepl("^s", colnames(testCounts))
    cObjSng <- spCounts(testCounts[,s], testErcc[,s])
    
    #setup expected data
    expected1 <- 1L
    expected2 <- 233L
    
    #run function
    output <- spTopMax(cObjSng, 250)
    
    #test
    expect_identical(expected1, output[1])
    expect_identical(expected2, output[length(output)])
    
})

##run test pearsonsDist
test_that("check that the pearsonsDist function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    s <- grepl("^s", colnames(testCounts))
    cObjSng <- spCounts(testCounts[,s], testErcc[,s])
    select <- 1:5

    #setup expected data
    expectHead <- c(0,1,1,1,1,0)
    expectTail <- rep(0, 6)
    
    #run function
    output <- pearsonsDist(cObjSng, select)
    
    #test
    expect_equivalent(expectHead, round(head(output)))
    expect_equivalent(expectTail, round(tail(output)))
    
})

##run test runTsne
test_that("check that the runTsne function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    my.dist <- as.dist(matrix(
        c(
            0, 9, 19,
            1, 10, 20,
            2, 11, 21),
            nrow=3
    ))
    
    dims <- 2
    theta <- 0
    initial_dims <- 3
    max_iter <- 5
    perplexity <- 0.1
    seed <- 11
    is_distance <- TRUE
    
    #setup expected data
    #no expected at the moment for this function
    
    #run function
    output <- runTsne(
        my.dist,
        dims,
        theta,
        initial_dims,
        max_iter,
        perplexity,
        seed,
        is_distance
    )
    
    #test
    expect_length(output, 6)
    
})

##run test runMclust
test_that("check that the runMclust function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    my.tsne <- matrix(
        c(
            0, 1,
            1, 2,
            4, 5,
            100, 111,
            101, 99,
            121, 102
        ),
        nrow=6,
        byrow=TRUE
    )
    
    Gmax <- 2
    seed <- 11
    
    #setup expected data
    exp_classification <- c(
        "A1", "A1", "A1",
        "B1", "B1", "B1"
    )
    exp_uncertainty <- rep(0, 6)
    
    #run function
    output <- runMclust(my.tsne, Gmax, seed)
    out_classes <- output[[1]]
    out_uncertainty <- output[[2]]
    
    #test
    expect_identical(out_classes, exp_classification)
    expect_identical(out_uncertainty, exp_uncertainty)
    
})

##run test averageGroupExpression
test_that("check that the averageGroupExpression function outputs the expected
result", {
    
    ###TEST1####
    #prepare normal input data
    s <- grepl("^s", colnames(testCounts))
    cObjSng <- spCounts(testCounts[,s], testErcc[,s])
    classes <- getData(testUns, "classification")

    #setup expected data
    expectedFirstRow <- c(5954, 998, 1651, 1054)
    expectedLastRow <- c(3603, 289, 1530, 753)

    names(expectedFirstRow) <- unique(classes)
    names(expectedLastRow) <- unique(classes)


    #run function
    output <- averageGroupExpression(cObjSng, classes, weighted=FALSE)
    
    #test
    expect_identical(expectedFirstRow, round(output[1,]))
    expect_identical(expectedLastRow, round(output[nrow(output),]))
    expect_false(any(is.na(output)))
})


##run test tsneGroupMeans
test_that("check that the tsneGroupMeans function outputs the expected result",{
    
    ###TEST1####
    #prepare normal input data
    classes <- getData(testUns, "classification")
    tsne <- getData(testUns, "tsne")
    
    #setup expected data
    expected <- data.frame(
        classification=c("A1", "B1", "C1", "D1"),
        x=c(-33,44,-15,4),
        y=c(-2,-6,-40,47)
    )
    
    #run function
    output <- tsneGroupMeans(tsne, classes)
    
    #round output
    output$x <- round(output$x)
    output$y <- round(output$y)
    
    #test
    expect_identical(expected, output)
    
})
