#context("spPlotSwarm")

s <- grepl("^s", colnames(testCounts))
cObjSng <- spCounts(testCounts[ ,s], testErcc[ ,s])
cObjMul <- spCounts(testCounts[ ,!s], testErcc[ ,!s])
uObj <- testUns
sObj <- testSwa

##run test plotSwarmGraph
test_that("check that plotSwarmGraph outputs the expected result", {
  
  #run function
  output <- plotSwarmGraph(sObj, uObj)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmBarBase
test_that("check that plotSwarmBarBase outputs the expected result", {
  
  #run function
  output <- plotSwarmBarBase(sObj)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmEdgeBar
test_that("check that plotSwarmEdgeBar outputs the expected result", {
  
  #run function
  output <- plotSwarmEdgeBar(sObj)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmPbar
test_that("check that plotSwarmPbar outputs the expected result", {
  
  #run function
  output <- plotSwarmPbar(sObj)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})

##run test plotSwarmHeat
test_that("check that plotSwarmHeat outputs the expected result", {
  
  #run function
  output <- plotSwarmHeat(sObj)
  
  #test
  expect_type(output, "list")
  expect_is(output, c("gg", "ggplot"))
})
