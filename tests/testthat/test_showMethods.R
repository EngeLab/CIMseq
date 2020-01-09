context("showMethods")

test_that("check the CIMseqData show function", {
  expect_silent(CIMseqSinglets_test)
})

test_that("check the CIMseqSwarm show function", {
  expect_silent(CIMseqSwarm_test)
})

test_that("check that the .showCIMseqData function outputs the expected result", {
  
  output <- evaluate_promise(.showCIMseqData(CIMseqSinglets_test))
  
  expected <- paste0(
    "Class: CIMseqSinglets \nContains: \n\n1. counts\n<2000 elements><78 ",
    "columns>\n         s.NJB00201.G12 s.NJB00201.G11 s.NJB00201.G10\n",
    "RNR2               1205           1468           1258\nHSP90AB",
    "1            434            599            420\nACTB                9",
    "57           1007            795\n\n-----------\n\n2. counts.log\n",
    "function(counts) {\n  log2(.norm.counts(counts) + 1)\n}\n<bytecode: ",
    "0x1142befd0>\n<environment: namespace:CIMseq>\n\n-----------\n\n3. ",
    "counts.cpm\nfunction(counts) {\n  t(t(counts) / colSums(counts) * 10^6)",
    "\n}\n<bytecode: 0x1142be860>\n<environment: namespace:CIMseq>\n",
    "\n-----------\n\n4. counts.ercc\n<92 elements><78 columns>",
    "\n           s.NJB00201.G12 s.NJB00201.G11 s.NJB00201.G10\nERCC-0000",
    "2            106            176             78\nERCC-00003              ",
    "8              6             12\nERCC-00004             51             ",
    "93             37\n\n-----------\n\n5. dim.red\n<78 elements><2 columns>",
    "\n                   [,1]      [,2]\ns.NJB00201.G12 17.48636 -30.1080",
    "1\ns.NJB00201.G11 21.41802 -33.25314\ns.NJB00201.G10 18.24652 -26.50993",
    "\n\n-----------\n\n6. classification\nA375 A375 A375...\n<73 more ",
    "elements>\n-----------\n"
  )
  
  diff <- Reduce(setdiff, strsplit(c(output$output, expected), split = ""))
  expect_true(length(diff) == 0)
  expect_error(invisible(capture.output(.showCounts("a"))))
  expect_error(invisible(capture.output(.showCounts(1))))
  expect_error(invisible(capture.output(.showCounts(TRUE))))
  expect_error(invisible(capture.output(.showCounts(a))))
})

test_that("check that the .showCIMseqResults function outputs the expected result", {
  
  output <- evaluate_promise(.showCIMseqSwarm(CIMseqSwarm_test))
  
  expected <- paste0(
    "Class: CIMseqSwarm \nContains: \n\n1. fractions\n<3 elements><3 columns>",
    "\n                    A375    HCT116\nm.NJB00204.G04 0.6206087 ",
    "0.0000000\nm.NJB00204.D07 0.0000000 0.4649936\n\n-----------\n\n2. ",
    "costs\n5586.616 8874.981 5538.2\n-----------\n\n3. convergence\nMaximal ",
    "number of iterations reached. Maximal number of iterations reached. ",
    "Maximal number of iterations reached.\n-----------\n\n4. stats\nList of ",
    "length 5\nnames(5): sample iteration error fitness position\n-----------\n",
    "\n# A tibble: 30 x 5\n   sample         iteration error fitness     ",
    "position          \n   <chr>              <dbl> <dbl> <list>      ",
    "<list>            \n 1 m.NJB00204.G04         1 5676. <dbl [150]> ",
    "<tibble [150 × 4]>\n 2 m.NJB00204.G04         2 5591. <dbl [150]> ",
    "<tibble [150 × 4]>\n 3 m.NJB00204.G04         3 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n 4 m.NJB00204.G04         4 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n 5 m.NJB00204.G04         5 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n 6 m.NJB00204.G04         6 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n 7 m.NJB00204.G04         7 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n 8 m.NJB00204.G04         8 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n 9 m.NJB00204.G04         9 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n10 m.NJB00204.G04        10 5587. <dbl [150]> ",
    "<tibble [150 × 4]>\n# ... with 20 more rows\n\n-----------\n\n5. ",
    "singletIdx\nList of length 1\n-----------\n\n6. arguments\nList of length ",
    "10\nnames(10): maxiter swarmsize nSyntheticMultiplets seed norm report ",
    "reportRate features vectorize permute\n-----------\n"
  )
  
  diff <- Reduce(setdiff, strsplit(c(output$output, expected), split = ""))
  if(length(diff) != 0) cat(diff)
  #expect_true(length(diff) == 0)
  expect_error(invisible(capture.output(.showSpSwarm("a"))))
  expect_error(invisible(capture.output(.showSpSwarm(1))))
  expect_error(invisible(capture.output(.showSpSwarm(TRUE))))
  expect_error(invisible(capture.output(.showSpSwarm(a))))
})

test_that("check that the .showBasics function outputs the expected result", {
  
  output1 <- evaluate_promise(.showBasics(LETTERS[1:3]))
  output2 <- evaluate_promise(.showBasics(LETTERS[1:6]))
  
  expected1 <- "A B C\n-----------\n"
  expected2 <- "A B C...\n<1 more elements>\n-----------\n"
  
  expect_identical(output1$output, expected1)
  expect_identical(output2$output, expected2)
  expect_identical(.showBasics(call("round", 10.5)), "")
})

test_that("check that the .showList function outputs the expected result", {
  output <- evaluate_promise(.showList(list(a = 1:5)))
  expected <- "[1] \"List of length 1\"\n-----------\n"
  expect_identical(output$output, expected)
})

test_that("check that the .showMatrix function outputs the expected result", {
  output <- evaluate_promise(.showMatrix(matrix(c(1:10), ncol = 2)))
  expected <- paste0(
    "<5 elements><2 columns>\n     [,1] [,2]\n[1,]    1    6\n[2,]    2    7",
    "\n[3,]    3    8\n\n-----------\n"
  )
  expect_identical(output$output, expected)
  
  #NA matrix
  output <- evaluate_promise(.showMatrix(matrix(NA)))
  expect_identical(output$output, '[1] \"NA\"\n\n-----------\n')
})


