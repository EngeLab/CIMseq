context("showMethods")

test_that("check that the .showCounts function outputs the expected result", {
  
  output <- evaluate_promise(.showCounts(test_spCountsSng))
  
  expected <- paste0(
    "Class: spCounts\nContains: \n1. counts\n<2000 ",
    "elements><81 columns>\n-----------\n\n2. counts.log\nfunction(counts) ",
    "{\n  log2(.norm.counts(counts) + 1)\n}\n<bytecode: ",
    "0x11132fdd8>\n<environment: namespace:sp.scRNAseq>\n\n-----------\n\n3. ",
    "counts.cpm\nfunction(counts) {\n  t(t(counts) / colSums(counts) * ",
    "10^6)\n}\n<bytecode: 0x1113314a8>\n<environment: ",
    "namespace:sp.scRNAseq>\n\n-----------\n\n4. counts.ercc\n<92 ",
    "elements><81 columns>\n-----------\n"
  )
  
  diff <- Reduce(setdiff, strsplit(c(output$output, expected), split = ""))
  if(length(diff) != 0) cat(diff)
  #expect_true(length(diff) == 0)
  expect_error(invisible(capture.output(.showCounts("a"))))
  expect_error(invisible(capture.output(.showCounts(1))))
  expect_error(invisible(capture.output(.showCounts(TRUE))))
  expect_error(invisible(capture.output(.showCounts(a))))
})

test_that("check that the .showUnsupervised function outputs the expected result", {
  
  output <- evaluate_promise(.showUnsupervised(test_spUnsupervised))
  
  expected <- paste0(
    "Class: spUnsupervised\nContains: \n1. tsne\n<81 elements><2 ",
    "columns>\n-----------\n\n2. tsneMeans\n<3 elements><3 ",
    "columns>\n-----------\n\n3. groupMeans\n<2000 elements><3 ",
    "columns>\n-----------\n\n4. classification\nA375 A375 A375 A375 HCT116 ",
    "HCT116...\n<76 more elements>\n-----------\n\n5. uncertainty\n0 0 0 0 0 ",
    "0...\n<76 more elements>\n-----------\n\n6. selectInd\n1 3 63 4 8 ",
    "55...\n<1995 more elements>\n-----------\n"
  )
  
  diff <- Reduce(setdiff, strsplit(c(output$output, expected), split = ""))
  if(length(diff) != 0) cat(diff)
  expect_true(length(diff) == 0)
  expect_error(invisible(capture.output(.showUnsupervised("a"))))
  expect_error(invisible(capture.output(.showUnsupervised(1))))
  expect_error(invisible(capture.output(.showUnsupervised(TRUE))))
  expect_error(invisible(capture.output(.showUnsupervised(a))))
})

test_that("check that the .showSpSwarm function outputs the expected result", {
  
  output <- evaluate_promise(.showSpSwarm(test_spSwarm))
  
  expected <- paste0(
    "Class: spSwarm\nContains: \n1. spSwarm\n<3 elements><3 ",
    "columns>\n-----------\n\n2. costs\n5586.71 8163.148 ",
    "5532.521\n-----------\n\n3. convergence\nMaximal number of iterations ",
    "reached. Maximal number of iterations reached. Maximal number of ",
    "iterations reached.\n-----------\n\n4. stats\n<30 elements><5 ",
    "columns>\n-----------\n\n5. singletIdx\nList of length ",
    "400\n-----------\n\n6. arguments\nList of length 10\nnames(10): maxiter ",
    "swarmsize nSyntheticMultiplets seed norm report reportRate selectInd ",
    "vectorize permute\n-----------\n"
  )
  
  diff <- Reduce(setdiff, strsplit(c(output$output, expected), split = ""))
  if(length(diff) != 0) cat(diff)
  expect_true(length(diff) == 0)
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
  expected <- "List of length 1\nnames(1): a\n-----------\n"
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


