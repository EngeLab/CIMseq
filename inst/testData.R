library(sp.scRNAseqTesting)

#expData

load('inst/counts.rda')
load('inst/counts.ercc.rda')

expData <- spCounts(counts, counts.ercc, '1000102901')
save(expData, file='data/expData.rda', compress='bzip2')

#unit test and vignettes data
load('data/syntheticData.rda')

#minimize cells
s.A1 <- syntheticData[ ,grepl("s.A1", colnames(syntheticData))][ ,1:80]
s.B1 <- syntheticData[ ,grepl("s.B1", colnames(syntheticData))][ ,1:80]
s.C1 <- syntheticData[ ,grepl("s.C1", colnames(syntheticData))][ ,1:80]
s.D1 <- syntheticData[ ,grepl("s.D1", colnames(syntheticData))][ ,1:80]
s.E1 <- syntheticData[ ,grepl("s.E1", colnames(syntheticData))][ ,1:80]
s.F1 <- syntheticData[ ,grepl("s.F1", colnames(syntheticData))][ ,1:80]
s.G1 <- syntheticData[ ,grepl("s.G1", colnames(syntheticData))][ ,1:80]
s.H1 <- syntheticData[ ,grepl("s.H1", colnames(syntheticData))][ ,1:80]
s.I1 <- syntheticData[ ,grepl("s.I1", colnames(syntheticData))][ ,1:80]
s.J1 <- syntheticData[ ,grepl("s.J1", colnames(syntheticData))][ ,1:80]

#add multuplets
m.A1B1C1D1 <- syntheticData[ ,'m.A1B1C1D1']
m.G1H1I1J1 <- syntheticData[ ,'m.G1H1I1J1']
m.A1B1 <- syntheticData[ ,'m.A1B1']

counts <- cbind(s.A1, s.B1, s.C1, s.D1, s.E1, s.F1, s.G1, s.H1, s.I1, s.J1, m.A1B1, m.A1B1C1D1, m.G1H1I1J1)

#minimise genes
.ntopMax <- function(data, n) {
    rv = apply(data, 1, max)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
}

select <- .ntopMax(counts, 250)
testData <- counts[select, ]
save(testData, file="data/testData.rda", compress="bzip2")
