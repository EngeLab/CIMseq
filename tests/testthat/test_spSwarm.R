#context("spSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])
uObj <- testUns
sObj <- testSwa

##run test .defineImport
test_that("check that the . function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    
    #setup expected data
    
    #run function

    #test

})

