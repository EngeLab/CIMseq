#context("spUnsupervised")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])

##run test spTopVar
test_that("check that spTopVar outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    expected1 <- 191L
    expected2 <- 25L
    
    #run function
    output <- spTopVar(cObjSng, 10)
    
    #test
    expect_identical(expected1, output[1])
    expect_identical(expected2, output[length(output)])
})

##run test spTopMax
test_that("check that spTopMax outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    expected1 <- 1L
    expected2 <- 107L
    
    #run function
    output <- spTopMax(cObjSng, 10)
    
    #test
    expect_identical(expected1, output[1])
    expect_identical(expected2, output[length(output)])
})

##run test pearsonsDist
test_that("check that pearsonsDist outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    select <- 1:5

    #setup expected data
    expectHead <- c(1, 2, 0, 1, 0, 2)
    expectTail <- c(0, 1, 0, 0, 0, 0)
    
    #run function
    output <- pearsonsDist(cObjSng, select)
    
    #test
    expect_equivalent(expectHead, round(head(output)))
    expect_equivalent(expectTail, round(tail(output)))
})

##run test runTsne
test_that("check that runTsne outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    my.dist <- as.dist(matrix(
        c(
            0, 9, 19,
            1, 10, 20,
            2, 11, 21),
            nrow = 3
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
test_that("check that runMclust outputs the expected result", {
    
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
test_that("check that averageGroupExpression outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    classes <- getData(testUns, "classification")

    #setup expected data
    expectedFirstRow <- c(4624, 598, 1505, 1920)
    expectedLastRow <- c(3921, 827, 3327, 6222)

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
test_that("check that tsneGroupMeans outputs the expected result",{
    
    ###TEST1####
    #prepare normal input data
    classes <- getData(testUns, "classification")
    tsne <- getData(testUns, "tsne")
    
    #setup expected data
    expected1 <- c(-42, 6)
    expected2 <- c(-9, -43)
    
    #run function
    output <- tsneGroupMeans(tsne, classes)
    
    #test
    expect_identical(expected1, as.numeric(round(output[1, 2:3])))
    expect_identical(expected2, as.numeric(round(output[nrow(output), 2:3])))
})

##run test erccPerClass
test_that("check that erccPerClass outputs the expected result",{
    
    ###TEST1####
    #setup expected data
    expected1 <- tibble::tibble(
        class = c("A1", "B1", "C1", "D1"),
        medianFracErcc = rep(1, 4)
    )
    
    #run function
    output <- erccPerClass(cObjSng, cObjMul, testUns)
    output$medianFracErcc <- round(output$medianFracErcc)
    
    #test
    expect_identical(expected1, output)
    expect_type(output$class, "character")
    expect_type(output$medianFracErcc, "double")
    expect_false(any(is.na(output$medianFracErcc)))
})

