#run from package root with: source('inst/rawData/testData.R')
packages <- c("CIMseq", "tidyverse", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

load('data/testMeta.rda')
load('data/testCounts.rda')

s <- grepl("^s", colnames(testCounts))
ercc <- grepl("^ERCC\\-[0-9]*$", rownames(testCounts))
singlets <- testCounts[!ercc, s]
singletsERCC <- testCounts[ercc, s]

#feature selection and dim red.
select <- selectTopMax(singlets, 2000)
pdist <- pearsonsDist(singlets, select)
tsne <- runTsne(pdist)
cObjSng <- CIMseqSinglets(
  singlets, singletsERCC, tsne, 
  testMeta$cellTypes[match(colnames(singlets), testMeta$sample)]
)
cObjMul <- CIMseqMultiplets(testCounts[!ercc, !s], testCounts[ercc, !s], select)

##spSwarm
future::plan(multiprocess)
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 10, swarmsize = 150,
  nSyntheticMultiplets = 400, report = TRUE, reportRate = 1
)

#rename
CIMseqSinglets_test <- cObjSng
CIMseqMultiplets_test <- cObjMul
CIMseqSwarm_test <- sObj
save(
  CIMseqSinglets_test, CIMseqMultiplets_test, CIMseqSwarm_test,
  file = "data/testData.rda", compress = "bzip2"
)
