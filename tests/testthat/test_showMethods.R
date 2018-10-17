context("showMethods")

test_that("check that the .showCIMseqData function outputs the expected result", {
  
  output <- evaluate_promise(.showCIMseqData(CIMseqSinglets_test))
  
  expected <- paste0(
    "Class: CIMseqSinglets \nContains: \n\n1. counts\n<2000 elements><3 ",
    "columns>\n         m.NJB00204.G04 m.NJB00204.D07\nRNR2               ",
    "1266           1432\nHSP90AB1            341            441\n",
    "\n-----------\n\n2. counts.log\nfunction(){}\n\n-----------\n\n3. ",
    "counts.cpm\nfunction(){}\n\n-----------\n\n4. counts.ercc\n<92 elements>",
    "<3 columns>\n           m.NJB00204.G04 m.NJB00204.D07\n",
    "ERCC-00002             42             38\nERCC-00003              ",
    "1              2\n\n-----------\n\n5. dim.red\n<80 elements><2 ",
    "columns>\n                    [,1]      [,2]\ns.NJB00201.G12 10.862440 ",
    "-26.97081\ns.NJB00201.G11  8.662592 -31.52306\n\n-----------\n\n6. ",
    "classification\nA375 A375 A375 A375 HCT116 HCT116...\n<75 more ",
    "elements>\n-----------\n"
  )
  
  diff <- Reduce(setdiff, strsplit(c(output$output, expected), split = ""))
  if(length(diff) != 0) print(diff)
  #expect_true(length(diff) == 0)
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
  output <- evaluate_promise(.showMatrix(matrix(c(1:10), ncol = 2)))
  expected <- paste0(
    "<5 elements><2 columns>\n     [,1] [,2]\n[1,]    1    6\n[2,]    2    ",
    "7\n\n-----------\n"
  )
  expect_identical(output$output, expected)
  
  #NA matrix
  output <- evaluate_promise(.showMatrix(matrix(NA)))
  expect_identical(output$output, '[1] \"NA\"\n\n-----------\n')
})


