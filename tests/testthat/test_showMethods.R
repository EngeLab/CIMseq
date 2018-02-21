#context("showMethods")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[, s], testErcc[, s])
cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
uObj <- testUns
sObj <- testSwa

test_that("check that the .showCounts function outputs the expected result", {
  
  output <- evaluate_promise(.showCounts(cObjSng))
  
  expected <- paste0(
    "Class: spCounts\nContains: \n1. counts\n<250 elements><340 columns>\n",
    "-----------\n\n2. counts.log\n<250 elements><340 columns>\n-----------\n",
    "\n3. counts.cpm\n<250 elements><340 columns>\n-----------\n\n4. ",
    "counts.ercc\n<2 elements><340 columns>\n-----------\n"
  )
  
  expect_identical(output$output, expected)
  expect_error(invisible(capture.output(.showCounts("a"))))
  expect_error(invisible(capture.output(.showCounts(1))))
  expect_error(invisible(capture.output(.showCounts(TRUE))))
  expect_error(invisible(capture.output(.showCounts(a))))
})

test_that("check that the .showUnsupervised function outputs the expected result", {
  
  output <- evaluate_promise(.showUnsupervised(uObj))
  
  expected <- paste0(
    "Class: spUnsupervised\nContains: \n1. tsne\n<340 elements><2 columns>\n",
    "-----------\n\n2. tsneMeans\n<4 elements><3 columns>\n-----------\n\n3. ",
    "groupMeans\n<250 elements><4 columns>\n-----------\n\n4. classification\n",
    "A1 A1 A1 A1 A1 A1...\n<335 more elements>\n-----------\n\n5. uncertainty",
    "\n0 0 0 0 0 0...\n<335 more elements>\n-----------\n\n6. selectInd\n1 23 ",
    "52 24 42 38...\n<245 more elements>\n-----------\n"
  )
  
  expect_identical(output$output, expected)
  expect_error(invisible(capture.output(.showUnsupervised("a"))))
  expect_error(invisible(capture.output(.showUnsupervised(1))))
  expect_error(invisible(capture.output(.showUnsupervised(TRUE))))
  expect_error(invisible(capture.output(.showUnsupervised(a))))
})

test_that("check that the .showSpSwarm function outputs the expected result", {
  
  output <- evaluate_promise(.showSpSwarm(sObj))
  
  expected <- paste0(
    "Class: spSwarm\nContains: \n1. spSwarm\n<2 elements><4 columns>\n",
    "-----------\n\n2. costs\n492647.2 451206.8\n-----------\n\n3. convergence",
    "\nMaximal number of iterations reached. Maximal number of iterations ",
    "reached.\n-----------\n\n4. stats\nList of length 0\n5. arguments\nList ",
    "of length 2\nnames(2): maxiter swarmsize"
  )
  
  expect_identical(output$output, expected)
  expect_error(invisible(capture.output(.showSpSwarm("a"))))
  expect_error(invisible(capture.output(.showSpSwarm(1))))
  expect_error(invisible(capture.output(.showSpSwarm(TRUE))))
  expect_error(invisible(capture.output(.showSpSwarm(a))))
})

test_that("check that the .showBasics function outputs the expected result", {
  
  output1 <- evaluate_promise(.showBasics(LETTERS[1:3]))
  output2 <- evaluate_promise(.showBasics(LETTERS[1:6]))
  
  expected1 <- "A B C\n-----------\n"
  expected2 <- "A B C D E F...\n<1 more elements>\n-----------\n"
  
  expect_identical(output1$output, expected1)
  expect_identical(output2$output, expected2)
  expect_identical(.showBasics(call("round", 10.5)), "")
})

test_that("check that the .showList function outputs the expected result", {
  output <- evaluate_promise(.showList(list(a = 1:5)))
  expected <- "List of length 1\nnames(1): a"
  expect_identical(output$output, expected)
})

test_that("check that the .showMatrix function outputs the expected result", {
  output <- evaluate_promise(.showMatrix(matrix(c(1:5), ncol = 1)))
  expected <- "<5 elements><1 columns>\n-----------\n"
  expect_identical(output$output, expected)
  
  #NA matrix
  output <- evaluate_promise(.showMatrix(matrix(NA)))
  expect_identical(output$output, '[1] \"NA\"\n\n-----------\n')
})


