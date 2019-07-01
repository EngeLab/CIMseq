context("plotSwarm")

##run test plotSwarmGraph
test_that("check that plotSwarmGraph outputs the expected result", {
  
  #run function
  output <- plotSwarmGraph(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmBarBase
test_that("check that plotSwarmBarBase outputs the expected result", {
  
  #run function
  output <- plotSwarmBarBase(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmEdgeBar
test_that("check that plotSwarmEdgeBar outputs the expected result", {
  
  #run function
  output <- plotSwarmEdgeBar(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmPbar
test_that("check that plotSwarmPbar outputs the expected result", {
  
  #run function
  output <- plotSwarmPbar(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmHeat
test_that("check that plotSwarmHeat outputs the expected result", {
  
  #run function
  output <- plotSwarmHeat(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
  )
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
