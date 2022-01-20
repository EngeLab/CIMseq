## ----style, echo = FALSE, results = 'asis', eval = TRUE-----------------------
BiocStyle::markdown()
packages <- c("CIMseq", "printr", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

## -----------------------------------------------------------------------------
counts.sng[1:2, 1:2] #example
dim(counts.sng) #genes and sample numbers
dim(counts.mul) #genes and sample numbers
class(counts.sng) #should be a matrix

## ---- eval = TRUE-------------------------------------------------------------
ercc.s <- grepl("^ERCC\\-[0-9]*$", rownames(counts.sng))
singlets <- counts.sng[!ercc.s, ]
singletERCC <- counts.sng[ercc.s, ]
ercc.m <- grepl("^ERCC\\-[0-9]*$", rownames(counts.mul))
multiplets <- counts.mul[!ercc.m, ]
multipletERCC <- counts.mul[ercc.m, ]

## ---- eval=TRUE---------------------------------------------------------------
#extract pre-calculated data
dim.red <- getData(CIMseqSinglets_test, "dim.red")
classes <- getData(CIMseqSinglets_test, "classification")
selected <- getData(CIMseqMultiplets_test, "features")

#build CIMseq objects
cObjSng <- CIMseqSinglets(counts=singlets, counts.ercc=singletERCC, classification=classes, dim.red=dim.red)
cObjMul <- CIMseqMultiplets(counts=multiplets, counts.ercc=multipletERCC, features=selected)

## ---- eval = TRUE-------------------------------------------------------------
counts <- getData(cObjSng, "counts")

## -----------------------------------------------------------------------------
empty <- CIMseqSinglets()
getData(empty, "counts") <- getData(cObjSng, "counts")

## ---- fig.align='center', fig.height=8, fig.width=10, eval=TRUE, message=FALSE----
plotCountsERCC(cObjSng, cObjMul)

## ---- fig.align='center', fig.height=8, fig.width=10, eval=TRUE, message=FALSE----
plotCountsMarkers(cObjSng, cObjMul, markers = c("CD74", "ANXA3"))

## ---- eval = FALSE------------------------------------------------------------
#  # For single-threaded, use "plan(sequential)" instead
#  plan(multicore)
#  sObj <- CIMseqSwarm(
#    cObjSng, cObjMul, swarmsize = 50,
#    nSyntheticMultiplets = 200, maxiter = 10, seed=123
#  )

## ---- echo = FALSE------------------------------------------------------------
sObj <- CIMseqSwarm_test

## ---- eval = TRUE-------------------------------------------------------------
calculateEdgeStats(sObj, cObjSng, cObjMul)

## ---- fig.align='center', fig.height=8, fig.width=10, eval=TRUE, message=FALSE----
plotSwarmCircos(sObj, cObjSng, cObjMul, weightCut=10, maxCellsPerMultiplet=4, alpha=Inf, h.ratio=0.9, depleted=F)

## ---- eval = TRUE-------------------------------------------------------------
mult.pred <- adjustFractions(singlets=cObjSng, multiplets=cObjMul, swarm=sObj, binary=T, maxCellsPerMultiplet=4)
mult.truth <- as.matrix(annot[rownames(mult.pred), c('A375', 'HCT116', 'HOS')])
pred.table <- mult.pred == as.logical(mult.truth)
mean(pred.table)

barplot(100 * (1-sapply(split(pred.table, rowSums(mult.truth)), mean)), ylab='% Error', xlab='Multiplet complexity', main='Deconvolution error rates')

## -----------------------------------------------------------------------------
sessionInfo()

