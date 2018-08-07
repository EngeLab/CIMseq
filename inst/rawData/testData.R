#run from package root with: source('inst/rawData/testData.R')
packages <- c("sp.scRNAseq", "tidyverse", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

meta <- read_tsv('./inst/rawData/testMeta.txt')
counts <- read_tsv('./inst/rawData/testCounts.txt') %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  as.matrix()

s <- grepl("^s", colnames(counts))
ercc <- grepl("^ERCC\\-[0-9]*$", rownames(counts))

cObjSng <- spCounts(counts[!ercc, s], counts[ercc, s])
cObjMul <- spCounts(counts[!ercc, !s], counts[ercc, !s])

##spUnsupervised
uObj <- spUnsupervised(cObjSng)

#rename classes
midx <- match(rownames(getData(uObj, "tsne")), meta$sample)
classification(uObj) <- meta$cellTypes[midx]
groupMeans(uObj) <- averageGroupExpression(
  cObjSng, getData(uObj, "classification"), FALSE
)
tsneMeans(uObj) <- tsneGroupMeans(
  getData(uObj, "tsne"), getData(uObj, "classification")
)

##spSwarm
future::plan(multiprocess)
sObj <- spSwarm(
  cObjSng, cObjMul, uObj, maxiter = 10, swarmsize = 150,
  nSyntheticMultiplets = 400, saveSingletData = TRUE, report = TRUE,
  reportRate = 1
)

#rename
test_spCountsSng <- cObjSng
test_spCountsMul <- cObjMul
test_spUnsupervised <- uObj
test_spSwarm <- sObj
save(
  test_spCountsSng, test_spCountsMul, test_spUnsupervised, test_spSwarm,
  file = "data/testData.rda"
)
