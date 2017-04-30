#context("spSwarm")

##run test .defineImport
test_that("check that the .defineImport function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    cellTypes <- data.frame(
        X1 = 1,
        X2 = 10,
        X3 = 10^6,
        index = 1
    )
    
    slice <- data.frame(
        oneTwo = rowMeans(cellTypes[ ,1:2]),
        index=1
    )
    
    fractions <- rep(1/(ncol(cellTypes)-1), ncol(cellTypes)-1)
    
    #setup expected data
    cmd1 <- 'import pandas as pd'
    cmd2 <- 'import numpy as np'
    cmd3 <- 'import math'
    cmd4 <- 'cellTypes = pd.DataFrame(cellTypesDictionary)'
    cmd5 <- 'slice = pd.DataFrame(sliceDictionary)'
    cmd6 <- 'fractions = np.asarray(fractions)'
    cmd7 <- 'cellTypes.sort_values(by=\'index\', ascending=\'True\')'
    cmd8 <- 'slice.sort_values(by=\'index\', ascending=\'True\')'
    cmd9 <- 'cellTypes = cellTypes.drop(\'index\', 1)'
    cmd10 <- 'slice = slice.drop(\'index\', 1)'
    
    expected_cmds <- list(
        cmd1,
        cmd2,
        cmd3,
        cmd4,
        cmd5,
        cmd6,
        cmd7,
        cmd8,
        cmd9,
        cmd10
    )
    
    expected_X1 <- 1
    expected_X2 <- 10
    expected_X3 <- 10^6
    expected_oneTwo <- 5.5
    expected_fractions <- rep(round(1/3, digits=5), 3)
    
    #run function
    output <- .defineImport(cellTypes, slice, fractions)

    #test
    expect_identical(expected_cmds, output)
    expect_identical(expected_X1, python.get('cellTypes[\'X1\'].values.tolist()'))
    expect_identical(expected_X2, python.get('cellTypes[\'X2\'].values.tolist()'))
    expect_identical(expected_X3, python.get('cellTypes[\'X3\'].values.tolist()'))
    expect_identical(expected_oneTwo, python.get('slice[\'oneTwo\'].values.tolist()'))
    expect_equivalent(expected_fractions, python.get('fractions.tolist()'))

})

##run test .defineCost
test_that("check that the .defineCost function outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    cmd1 <- 'def makeSyntheticSlice(cellTypes, fractions):
        func = lambda x: sum(np.asarray(x) * np.asarray(fractions))
        return cellTypes.apply(func, axis=1)'
    
    cmd2 <- 'def distToSlice(fractions, *args):
        cellTypes, slice, col = args
        normFractions = fractions / sum(fractions)
        a = makeSyntheticSlice(cellTypes, normFractions)
        for i in a:
            if math.isnan(i):
                print "NaN in make synthetic slice!"
    
        diff = []
        for index in range(cellTypes.shape[0]):
            d = a.iloc[index] - slice.iloc[index, col]
            diff.append(abs(d))
        cost = sum(diff)
        return cost'
    
    expected <- list(cmd1, cmd2)
    
    #run function
    output <- .defineCost()
    
    #test
    expect_identical(expected, output)
    
})

##run test .con1
test_that("check that the .con1 function outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    expected <- 'def con1(fractions, *args):
        if sum(fractions) > 0.1:
            return 0
        else:
            return -1'
            
    
    #run function
    output <- .con1()
    
    #test
    expect_identical(expected, output)
    
})

##run test .definePySwarm
test_that("check that the .definePySwarm function outputs the expected result", {
    
    ###TEST1####
    #setup expected data
    expected <- 'def optimize(cellTypes, slice, fractions, col, maxiter, swarmsize, minstep, minfunc):
        lb = np.asarray([0] * len(fractions))
        ub = np.asarray([1] * len(fractions))
        name = slice.columns.values[col]
        args = (cellTypes, slice, col)
        xopt, fopt = pso(
            distToSlice,
            lb,
            ub,
            args=args,
            f_ieqcons=con1,
            maxiter=maxiter,
            swarmsize=swarmsize,
            minstep=minstep,
            minfunc=minfunc
        )
        dictionary = dict(zip(list(cellTypes.columns.values), xopt.tolist()))
        return { \'xopt\':dictionary, \'fopt\':fopt, \'name\':name }'
    
    
    #run function
    output <- .definePySwarm()
    
    #test
    expect_identical(expected, output)
    
})

##run test .runPySwarm
test_that("check that the .runPySwarm function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    cellTypes <- data.frame(
        X1 = 1,
        X2 = 10,
        X3 = 10^6,
        index = 1
    )
    
    slice <- data.frame(
        oneTwo = rowMeans(cellTypes[ ,1:2]),
        index=1
    )
    
    fractions <- rep(1/(ncol(cellTypes)-1), ncol(cellTypes)-1)
    limit <- "none"
    maxiter <- 100
    swarmsize <- 10
    minstep <- 1e-16
    minfunc <- 1e-16
    cores <- 1
    
    #setup expected data
    exp_currentXopt <- 0
    exp_currentFoptX1 <- c(X1 = 1)
    exp_currentFoptX2 <- c(X2 = 1)
    exp_currentFoptX3 <- c(X3 = 0)
    
    #run function
    .defineImport(cellTypes, slice, fractions)
    .defineCost()
    .con1()
    .definePySwarm()

    output <- .runPySwarm(
        cellTypes,
        slice,
        fractions,
        limit,
        maxiter,
        swarmsize,
        minstep,
        minfunc,
        cores
    )
    
    #test
    expect_equivalent(exp_currentXopt, output[['currentFopt']])
    expect_true("X1" %in% names(
        output[['currentXopt']][order(output[['currentXopt']], decreasing =TRUE)][1:2])
    )
    expect_true("X2" %in% names(
    output[['currentXopt']][order(output[['currentXopt']], decreasing =TRUE)][1:2])
    )

})

##run test .processResults
test_that("check that the .processResults function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    result <- list(
        currentFopt = 0,
        currentXopt = c(X2=1, X3=2, X1=3),
        name = "oneThree",
        currentFopt = 1,
        currentXopt = c(X2=10, X3=20, X1=30),
        name = "oneTwo"
    )
    
    #setup expected data
    expected <- data.frame(
        sampleName = c("oneThree", "oneTwo"),
        fopt = c(0,1),
        X1 = rep(1/2, 2),
        X2 = rep(1/6, 2),
        X3 = rep(1/3, 2),
        stringsAsFactors = FALSE
    )
    
    #run function
    output <- .processResults(result)
    
    #test
    expect_identical(output, expected)
    
})

##run test .multiHOTencoding
test_that("check that the .multiHOTencoding function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    optResult <- data.frame(
        sampleName = c("oneThree", "oneTwo"),
        fopt = c(0,1),
        X1 = rep(1/2, 2),
        X2 = rep(1/6, 2),
        X3 = rep(1/3, 2),
        stringsAsFactors = FALSE
    )
    
    cutoff <- 0.4
    
    #setup expected data
    expected <- data.frame(
        fopt = c(0,1),
        sampleName = c("oneThree", "oneTwo"),
        X1 = rep(1/2, 2),
        X2 = rep(0, 2),
        X3 = rep(0, 2),
        stringsAsFactors = FALSE
    )
    
    #run function
    output <- .multiHOTencoding(optResult, matrix(), cutoff)
    
    #test
    expect_identical(output, expected)
    
})

