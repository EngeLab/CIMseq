context("spSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[,s], testErcc[,s])
cObjMul <- spCounts(testCounts[,!s], testErcc[,!s])
uObj <- testUns
sObj <- testSwa

#Function to check if all elements in a vector are identical
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

##run test syntheticMultipletsFromCounts
test_that("check that syntheticMultipletsFromCounts outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- getData(cObjSng, "counts.cpm")
  classes <- getData(uObj, "classification")
  fractions <- rep(1, length(unique(classes)))
  seed <- 923
  
  #setup expected data
  nr <- 250L
  nc <- 1L
  first <- 5459L
  last <- 18599L
  
  #run function
  output <- syntheticMultipletsFromCounts(singlets, classes, fractions, seed)
  
  #test
  expect_is(output, "matrix")
  expect_silent(syntheticMultipletsFromCounts(singlets, classes, fractions, seed))
  expect_identical(ncol(output), nc)
  expect_identical(nrow(output), nr)
  expect_identical(output[[1]], first)
  expect_identical(output[[250]], last)
  expect_error(syntheticMultipletsFromCounts(classes, fractions, seed))
  expect_error(syntheticMultipletsFromCounts(singlets, fractions, seed))
  expect_error(syntheticMultipletsFromCounts(singlets, classes, seed))
  expect_error(syntheticMultipletsFromCounts(singlets, classes, fractions))
})

##run test costCalculation
test_that("check that costCalculation outputs the expected result", {
  
    ###TEST1####
    #prepare normal input data
    oneMultiplet <- getData(cObjMul, "counts.cpm")[, 1]
    singlets <- getData(cObjSng, "counts.cpm")
    classes <- getData(uObj, "classification")
    seed <- 923
    fractions.wrong <- rep(1, length(unique(classes)))
    fractions.right <- c(0.5, 0.5, 0, 0)
    
    syntheticMultiplets.wrong <- generateSyntheticMultiplets(
      singlets, classes, fractions.wrong, seed, 100
    )
    
    syntheticMultiplets.right <- generateSyntheticMultiplets(
      singlets, classes, fractions.right, seed, 100
    )
    
    #setup expected data
    expected <- 1639.734
    
    #run function
    output.wrong <- costCalculation(oneMultiplet, syntheticMultiplets.wrong)
    output.right <- costCalculation(oneMultiplet, syntheticMultiplets.right)
    
    #test
    expect_is(output.wrong, "numeric")
    expect_is(output.right, "numeric")
    expect_true(length(output.right) == 1)
    expect_true(length(output.wrong) == 1)
    expect_silent(costCalculation(oneMultiplet, syntheticMultiplets.wrong))
    expect_error(costCalculation(syntheticMultiplets.wrong))
    expect_error(costCalculation(oneMultiplet))
    expect_true(output.right < output.wrong)
    expect_true(all.equal(output.right, expected, tolerance = 1e-3))
})

##run test generateSyntheticMultiplets
test_that("check that generateSyntheticMultiplets outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  singlets <- getData(cObjSng, "counts.cpm")
  classes <- getData(uObj, "classification")
  seed <- 923
  fractions <- c(0.5, 0.5, 0, 0)
  
  #setup expected data
  nc <- 100
  nr <- 250
  
  #run function
  output <- generateSyntheticMultiplets(
    singlets, classes, fractions, seed, 100
  )
  
  #test
  expect_is(output, "matrix")
  expect_silent(generateSyntheticMultiplets(
    singlets, classes, fractions, seed, 100
  ))
  expect_error(generateSyntheticMultiplets(
    classes, fractions, seed, 100
  ))
  expect_error(generateSyntheticMultiplets(
    singlets, fractions, seed, 100
  ))
  expect_error(generateSyntheticMultiplets(
    singlets, classes, seed, 100
  ))
  expect_error(generateSyntheticMultiplets(
    singlets, classes, fractions, 100
  ))
  expect_error(generateSyntheticMultiplets(
    singlets, classes, fractions, seed
  ))
  expect_true(ncol(output) == nc)
  expect_true(nrow(output) == nr)
})

##run test cost.fn
test_that("check that cost.fn outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  oneMultiplet <- getData(cObjMul, "counts.cpm")[, 1]
  singlets <- getData(cObjSng, "counts.cpm")
  classes <- getData(uObj, "classification")
  seed <- 923
  fractions <- c(0.5, 0.5, 0, 0)
  
  syntheticMultiplets <- generateSyntheticMultiplets(
    singlets, classes, fractions, seed, 100
  )
  
  #setup expected data
  expected <- 1639.734
  expected.sumFrac0 <- 999999999
  
  #run function
  output.1 <- cost.fn(fractions, oneMultiplet, singlets, classes, seed, 100)
  output.2 <- cost.fn(rep(0, length(unique(classes))), oneMultiplet, singlets, classes, seed, 100)
  
  #test
  expect_true(all.equal(output.1, expected, tolerance = 1e-3))
  expect_identical(expected.sumFrac0, output.2)
})

##run test getMultipletsForEdge
test_that("check that getMultipletsForEdge outputs the expected result", {
  
    ###TEST1####
    #setup expected data
    #A1 and B1 should have an edge
    #I1 and J1 should have an edge
    expected1 <- tibble::tibble(
        multiplet = "m.A1B1.341",
        from = "A1",
        to = "B1"
    )
    expected2 <- tibble::tibble(
        multiplet = "m.C1D1.342",
        from = "C1",
        to = "D1"
    )
    expected3 <- tibble::tibble(
        multiplet = c("m.A1B1.341", "m.C1D1.342"),
        from = c("A1", "C1"),
        to = c("B1", "D1")
    )
    
    #run function
    output1 <- getMultipletsForEdge(sObj, 1/4, data.frame("A1", "B1"))
    output2 <- getMultipletsForEdge(sObj, 1/4, data.frame("C1", "D1"))
    output3 <- getMultipletsForEdge(
      sObj, 1/4, data.frame(c("A1", "C1"), c("B1", "D1"))
    )
    
    #test
    expect_identical(output1, expected1)
    expect_identical(output2, expected2)
    expect_identical(output3, expected3)
    
    ###TEST2####
    #prepare normal input data
    tmp <- sObj
    tmp@spSwarm <-
        data.frame(
            A1 = c(0.3, 0.3),
            B1 = c(0.3, 0.3),
            C1 = c(0.3, 0),
            D1 = c(0, 0.3),
            row.names=c("m.A1B1C1a", "m.A1B1C1b")
        )
    
    #setup expected data
    expected <- tibble::tibble(
        multiplet = rep(c("m.A1B1C1a", "m.A1B1C1b"), 3),
        from = c(rep("A1", 4), rep("B1", 2)),
        to = c(rep("B1", 2), rep(c("C1", "D1"), 2))
    )
    
    
    #run function
    output <- getMultipletsForEdge(
        tmp,
        1/4,
        data.frame(
            c("A1", "A1", "A1", "B1", "B1", "C1"),
            c("B1", "C1", "D1", "C1", "D1", "D1")
        )
    )
    
    #test
    expect_identical(output, expected)
})

##run test getEdgesForMultiplet
test_that("check that getEdgesForMultiplet outputs the expected result", {
  
    ###TEST1####
    #setup expected data
    expected1 <- tibble::tibble(
        multiplet = "m.A1B1",
        from = "A1",
        to = "B1"
    )
    expected2 <- tibble::tibble(
        multiplet = "m.C1D1",
        from = "C1",
        to = "D1"
    )

    #run function
    output1 <- getEdgesForMultiplet(sObj, 1/4, 'm.A1B1')
    output2 <- getEdgesForMultiplet(sObj, 1/4, 'm.C1D1')

    #test
    expect_identical(output1, expected1)
    expect_identical(output2, expected2)
})
