#context("spCounts")

##run test
test_that("check that the .norm.log.counts function outputs the expected result", {

    ###TEST1####
    #prepare normal input data
    input <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    
    #setup expected data
    expected <- matrix(
        c(
            0,
            log2(1000001),
            log2(1000001),
            0
        ),
        nrow=2,
        ncol=2
    )

    #run function
    output <- .norm.log.counts(input)

    #test
    expect_true(all.equal(expected, output))

})

test_that("check that the .sampleType function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    sampleType <- 'm.'
    counts <- matrix(
        c(0,1,1,0),
        nrow=2,
        ncol=2,
        dimnames=list(
            c("a", "b"),
            c("m.1", "s.1")
        )
    )
    
    #setup expected data
    expected <- c("Multuplet", "Singlet")
    
    #run function
    output <- .sampleType(sampleType, counts)
    
    #test
    expect_true(all.equal(expected, output))
    
})
