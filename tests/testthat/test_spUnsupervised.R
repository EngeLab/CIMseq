#context("spUnsupervised")

##run test tsneGroupMeans
test_that("check that the tsneGroupMeans function outputs the expected result",{
    
    ###TEST1####
    #prepare normal input data
    data <- matrix(
        c(
            1,2,1,2,
            2,3,2,3,
        ),
        nrow=4,
        ncol=2,
        byrow=TRUE
    )
    
    classes <- sort(rep(LETTERS[1:2], 2))
    
    #setup expected data
    expected <- data.frame(
        classification=c("A", "B"),
        x=c(1,2),
        y=c(2,3)
    )
    
    #run function
    output <- tsneGroupMeans(data, classes)
    
    #test
    expect_identical(expected, output)
    
})


##run test .ntopMax
test_that("check that the .ntopMax function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    data <- matrix(
        c(
            1,1,4,1,
            2,2,3,3,
            3,3,4,5,
            0,0,0,0
        ),
        nrow=4,
        ncol=4,
        byrow=TRUE
    )
    
    #setup expected data
    expected <- as.integer(c(3, 1))
    
    #run function
    output <- .ntopMax(data, 2)
    
    #test
    expect_identical(expected, output)
    
})


##run test averageGroupExpression
test_that("check that the averageGroupExpression function outputs the expected 
result", {
    
    ###TEST1####
    #prepare normal input data
    sng <- matrix(
        c(
            1,10,2,20,
            2,20,3,30,
            3,30,0,0,
            0,0,0,0
        ),
        nrow=4,
        ncol=4,
        byrow=TRUE
    )
    
    classes <- rep(LETTERS[1:2], 2)
    colnames(sng) <- classes

    #setup expected data
    expected <- matrix(
        c(
            1.5, 15,
            2.5, 25,
            1.5, 15,
            0, 0
        ),
        byrow=TRUE,
        nrow=4,
        dimnames=list(
            c(),
            c(unique(classes))
        )
    )
    
    #run function
    output <- averageGroupExpression(cObjSng, classes, weighted=FALSE)
    
    #test
    expect_identical(expected, output)
    expect_false(any(is.na(output)))
})

##run test .tsneGroupMeans
test_that("check that the .tsneGroupMeans function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    x <- matrix(
        c(
            1,2,2,3,3,0,0,0,
            2,1,3,2,0,3,0,0
        ),
        nrow=8,
        ncol=2,
    )
    
    class <- rep(LETTERS[1:2], 4)
    
    #setup expected data
    expected <- data.frame(
        classification = c(LETTERS[1:2]),
        meanX = c(1.5, 1.25),
        meanY = c(1.25, 1.5)
    )
    
    #run function
    output <- .tsneGroupMeans(x, class)
    
    #test
    expect_identical(expected, output)
    
})

##run test .distFunc
test_that("check that the .distFunc function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    x <- matrix(
        c(
            1,2,3,2,
            2,3,4,3,
            3,4,5,0,
            0,0,0,0
        ),
        nrow=4,
        ncol=4,
        byrow=TRUE
    )
    
    select <- c(1,2,3)
    sampleType <- c(rep("Singlet", 3), "Multuplet")
    
    #setup expected data
    expected <- as.dist(matrix(
        rep(0,9),
        nrow=3
    ))
    
    #run function
    output <- .distFunc(x, select, sampleType)
    
    #test
    expect_equivalent(expected, output)
    
})

##run test .runTsne
test_that("check that the .runTsne function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    my.dist <- as.dist(matrix(
        c(
            0, 9, 19,
            1, 10, 20,
            2, 11, 21),
        nrow=3
    ))
    
    k <- 2
    plot.callback <- NULL
    initial_dims <- 3
    max_iter <- 5
    perplexity <- 0.1
    seed <- 11
    
    #setup expected data
    #no expected at the moment for this function
    
    #run function
    output <- .runTsne(my.dist, k, initial_dims, max_iter, perplexity, seed, theta=0)
    
    #test
    expect_length(output, 6)
    
})

##run test .runMclust
test_that("check that the .runMclust function outputs the expected result", {
    
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
    
    #run function
    output <- .runMclust(seed, my.tsne, G=Gmax)
    
    #test
    expect_identical(output$classification, exp_classification)
    
})
