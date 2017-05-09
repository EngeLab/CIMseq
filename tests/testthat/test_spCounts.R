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

test_that("check that the .norm.counts function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    input <- matrix(c(0,1,1,0), nrow=2, ncol=2)
    
    #setup expected data
    expected <- matrix(
        c(
            0,
            1000001,
            1000001,
            0
        ),
        nrow=2,
        ncol=2
    )
    
    #run function
    output <- .norm.counts(input)
    
    #test
    expect_true(all.equal(expected, output))


})

test_that("check that the .inputCheckCounts function outputs the expected result", {
    
    ###TEST1####
    #non-conforming counts and counts.ercc
    #setup input
    counts <- matrix(1:10, ncol=10)
    counts.ercc <- matrix(1:11, ncol=11)
    
    #test
    expect_message(.inputCheckCounts(counts, counts.ercc))
    
    ###TEST2####
    #NA's present
    #setup input
    counts <- matrix(NA, ncol=10)
    counts.ercc <- matrix(NA, ncol=11)
    
    #test
    expect_message(.inputCheckCounts(counts, counts.ercc))
    
    
})
