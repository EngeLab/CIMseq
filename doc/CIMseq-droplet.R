## ----style, echo = FALSE, results = 'asis', eval = TRUE-----------------------
BiocStyle::markdown()
packages <- c("CIMseq", "printr", "future", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

## -----------------------------------------------------------------------------
SmallIntest.singlets.counts[1:2, 1:2] #example
dim(SmallIntest.singlets.counts) #genes and sample numbers
dim(SmallIntest.multiplets.counts) #genes and sample numbers
class(SmallIntest.singlets.counts) #should be a matrix

## ---- eval=TRUE---------------------------------------------------------------
#build CIMseq objects
cObjSng.si <- CIMseqSinglets(counts=SmallIntest.singlets.counts, classification=SmallIntest.singlets.classes, norm.to=10000)
cObjMul.si <- CIMseqMultiplets(counts=SmallIntest.multiplets.counts, features=SmallIntest.marker.genes, norm.to=10000)

## ---- eval = FALSE------------------------------------------------------------
#  # For single-threaded, use "plan(sequential)" instead
#  plan(multicore)
#  options(future.globals.maxSize= 1000000000)
#  sObj.si <- CIMseqSwarm(cObjSng.si, cObjMul.si, maxiter=100, swarmsize=110, nSyntheticMultiplets=200, seed=123)

## ---- eval = TRUE-------------------------------------------------------------
si.edges <- calculateEdgeStats(sObj.si, cObjSng.si, cObjMul.si, multiplet.factor=2, maxCellsPerMultiplet=2)
si.edges %>% filter(pval < 1e-4 & weight > 10) %>% arrange(desc(score)) %>% head(n=10)

## ---- echo = FALSE------------------------------------------------------------
classOrder <- c("Paneth", "Stem", "Transit amplifying", "Progenitor early", "Progenitor late-1", "Progenitor late-2", "Enterocyte", "Goblet", "Enteroendocrine", "Tuft")

## ---- fig.align='center', fig.height=8, fig.width=10, eval=TRUE, message=FALSE----
plotSwarmCircos(sObj.si, cObjSng.si, cObjMul.si, weightCut=10,
			 maxCellsPerMultiplet=2, alpha=1E-4, h.ratio=0.95,
			 depleted=F, multiplet.factor=2, classOrder=classOrder)

## -----------------------------------------------------------------------------
sessionInfo()

