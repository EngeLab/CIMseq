#context("Optimization")


##run test
#test_that("check that the .dist.to.slice function outputs the expected result", {
    
    ####TEST1####
    ##prepare normal input data
    #fake seq data
#    set.seed(11)
#    x <- c(sample(seq(0,10^6, 1), size=2000), x1 + 1, x1 - 1)
#    noiseV <- sample(seq(1,length(x), 1), size=(length(x)*0.3), replace=FALSE)
#    x[noiseV] <- 0
#    m <- as.data.frame(matrix(x, ncol=3))
    
    #make multuplets
#    test <- data.frame(rowMeans(m[ ,1:2]), rowMeans(m[ ,2:3]), rowMeans(m))
    
    ##setup expected data
#    expected1 <- 0
#    expected2 <- 0
#    expected3 <- 0
    
    ##run function
#    result1 <- .dist.to.slice(c(0.5, 0.5, 0), m, test[ ,1])
#   result2 <- .dist.to.slice(c(0.5, 0.5, 0), m, test[ ,2])
#   result3 <- .dist.to.slice(rep(1/3, 3), m, test[ ,2])
    
    ##test
#    expect_true(all.equal(expected1, result1))
#    expect_true(expected2 < result2)
#    expect_true(all.equal(expected3, result3))

#})
