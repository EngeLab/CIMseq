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