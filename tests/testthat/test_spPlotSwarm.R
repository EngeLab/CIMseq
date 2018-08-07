context("spPlotSwarm")

##run test plotSwarmGraph
test_that("check that plotSwarmGraph outputs the expected result", {
  
  #run function
  output <- plotSwarmGraph(test_spSwarm, test_spUnsupervised)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmBarBase
test_that("check that plotSwarmBarBase outputs the expected result", {
  
  #run function
  output <- plotSwarmBarBase(test_spSwarm)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmEdgeBar
test_that("check that plotSwarmEdgeBar outputs the expected result", {
  
  #run function
  output <- plotSwarmEdgeBar(test_spSwarm)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmPbar
test_that("check that plotSwarmPbar outputs the expected result", {
  
  #run function
  output <- plotSwarmPbar(test_spSwarm)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmHeat
test_that("check that plotSwarmHeat outputs the expected result", {
  
  #run function
  output <- plotSwarmHeat(test_spSwarm)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
