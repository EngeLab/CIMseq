#from package directory run with source('./inst/rawData/testData.R')

library(sp.scRNAseqTesting)

########unit test and vignettes data
#minimize cells
s.A1 <- syntheticData[, grepl("s.A1", colnames(syntheticData))][, 1:85]
s.B1 <- syntheticData[, grepl("s.B1", colnames(syntheticData))][, 1:85]
s.I1 <- syntheticData[, grepl("s.I1", colnames(syntheticData))][, 1:85]
s.J1 <- syntheticData[, grepl("s.J1", colnames(syntheticData))][, 1:85]

#add multuplets
#note that the names A1B1 and C1D1 are kept although, in the synthetic data,
#these are actually B1A1 and J1I1. This is necessary because during the
#spUnsupervised phase the names of the cell types are automatically set using a
#combination of letters and numbers. Each cell type will be names with A1, B1,
#C1, etc. Therefore, since there are 4 cell types in the data, we know their
#names ahead of time and setting the multiplet names to something else only
#makes a correct result look incorrect.

m.A1B1 <- syntheticData[, 'm.A1B1']
m.C1D1 <- syntheticData[, 'm.I1J1']

#make counts
counts <- cbind(
  s.A1,
  s.B1,
  s.I1,
  s.J1,
  m.A1B1,
  m.C1D1
)

#minimise genes
.ntopMax <- function(data, n) {
  rv = apply(data, 1, max)
  select = order(rv, decreasing = TRUE)[1:n]
  return(select)
}


select <- .ntopMax(counts, 250)
testCounts <- counts[select, ]
rownames(testCounts) <- sort(
  paste(
    rep(
      letters,
      10
    ),
    1:11,
    sep = ""
  )[1:nrow(testCounts)]
)

#make testErcc
s <- grepl("^s", colnames(expCounts))
s2 <- grepl("^s", colnames(testCounts))
singletsE <- expErcc[c(1, 2), s]
singletsE <- singletsE[, sample(
  1:ncol(expErcc[, s]),
  size = length(s2[s2]),
  replace = TRUE
)]

set.seed(2342536)
idx <- sample(1:ncol(singletsE), size = ceiling(ncol(singletsE) * 0.8), replace = FALSE)
m1 <- rowMeans(singletsE[, idx]) / 9
set.seed(54254)
idx <- sample(1:ncol(singletsE), size = ceiling(ncol(singletsE) * 0.8), replace = FALSE)
m2 <- rowMeans(singletsE[, idx]) / 9

multipletsE <- matrix(c(m1, m2), ncol = 2)
colnames(multipletsE) <- c("m.A1B1", "m.C1D1")

testErcc <- cbind(singletsE, multipletsE)

#make test spUnsupervised and spSwarm
cObjSng <- spCounts(testCounts[, s2], testErcc[, s2])
cObjMul <- spCounts(testCounts[, !s2], testErcc[, !s2])
testUns <- spUnsupervised(cObjSng, max = 250, max_iter = 1000)
cn <- estimateCells(cObjSng, cObjMul)

testSwa <- spSwarm(
  cObjMul,
  testUns,
  distFun = "dtsnCellNum",
  maxiter = 100,
  swarmsize = 500,
  cores = 2,
  cellNumbers = cn,
  e = 0.0025
)

#save
save(
  testErcc,
  testCounts,
  testUns,
  testSwa,
  file = "data/testData.rda",
  compress = "bzip2"
)
