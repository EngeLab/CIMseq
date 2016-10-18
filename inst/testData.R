library(sp.scRNAseqTesting)

#expData

load('inst/counts.rda')
load('inst/counts.ercc.rda')

expData <- spCounts(counts, counts.ercc, '1000102901')
save(expData, file='data/expData.rda', compress='bzip2')

#unit test and vignettes data

testData <- syntheticData[,1:1000]

A1 <- testData[ ,grepl("s.A1", colnames(testData))][1:2000,1:100]
B1 <- testData[ ,grepl("s.B1", colnames(testData))][1:2000,1:100]
C1 <- testData[ ,grepl("s.C1", colnames(testData))][1:2000,1:100]
D1 <- testData[ ,grepl("s.D1", colnames(testData))][1:2000,1:100]
E1 <- testData[ ,grepl("s.E1", colnames(testData))][1:2000,1:100]
F1 <- testData[ ,grepl("s.F1", colnames(testData))][1:2000,1:15]
G1 <- testData[ ,grepl("s.G1", colnames(testData))][1:2000,1:15]
H1 <- testData[ ,grepl("s.H1", colnames(testData))][1:2000,1:15]
I1 <- testData[ ,grepl("s.I1", colnames(testData))][1:2000,1:15]
J1 <- testData[ ,grepl("s.J1", colnames(testData))][1:2000,1:15]

mult <- c("m.A1B1", "m.B1C1", "m.C1D1", "m.D1E1")
testData <- cbind(A1, B1, C1, D1, E1, syntheticData[1:nrow(A1) ,mult])

#F1, G1, H1, I1, J1, syntheticData[1:nrow(A1) ,mult]
##ercc
ercc <- matrix(c(rep(rep(9, 15), rep(4.5, 3), 3), 2), nrow=2)

##spCounts
cObj <- spCounts(testData, ercc, "m.")

##spUnsupervised
uObj <- spUnsupervised(cObj, max=2000, max_iter=2000, Gmax=5)
spPlot(uObj, type="clusters")
table(getData(uObj, "classification"))

#spSwarm
sObj <- spSwarm(cObj, uObj, cores=4, swarmsize=250, maxiter=10)


##sanity check
cObj <- spCounts(syntheticData, matrix(), "m.")
uObj <- spUnsupervised(cObj, max=2000, max_iter = 1000)
sObj <- spSwarm(uObj, limit=2, cores=2, swarmsize=250, maxiter=10)
