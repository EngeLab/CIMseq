library(sp.scRNAseqTesting)

#expData
load('inst/counts.rda')
load('inst/counts.ercc.rda')

#expData <- spCounts(counts, counts.ercc, '1000102901')
#save(expData, file='data/expData.rda', compress='bzip2')

expCounts <- counts
expErcc <- counts.ercc

suffix <- ifelse(grepl('1000102901', colnames(counts)), "s", "m")
colnames(expCounts) <- paste(suffix, colnames(counts), sep='.')
colnames(expErcc) <- paste(suffix, colnames(counts.ercc), sep='.')

save(expCounts, file='data/expCounts.rda', compress='bzip2')
save(expErcc, file='data/expErcc.rda', compress='bzip2')

########unit test and vignettes data
#minimize cells
s.A1 <- syntheticData[ ,grepl("s.A1", colnames(syntheticData))][ ,1:85]
s.B1 <- syntheticData[ ,grepl("s.B1", colnames(syntheticData))][ ,1:85]
s.I1 <- syntheticData[ ,grepl("s.I1", colnames(syntheticData))][ ,1:85]
s.J1 <- syntheticData[ ,grepl("s.J1", colnames(syntheticData))][ ,1:85]

#add multuplets
m.A1B1 <- syntheticData[ ,'m.A1B1']
m.C1D1 <- syntheticData[ ,'m.I1J1']

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
    select = order(rv, decreasing=TRUE)[1:n]
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
        sep=""
    )[1:nrow(testCounts)]
)

#make testErcc
s <- grepl("^s", colnames(expCounts))
s2 <- grepl("^s", colnames(testCounts))
set.seed(129)
singletsE <- expErcc[,s]
singletsE <- singletsE[,sample(
    1:ncol(expErcc[,s]),
    size=length(s2[s2 == TRUE]),
    replace=TRUE
)]

multipletsE <-expErcc[,!s]
multipletsE <- multipletsE[,
    sample(1:ncol(expErcc[,!s]),
    size=length(s2[s2 == FALSE]),
    replace=TRUE
)]

testErcc <- matrix(c(singletsE, multipletsE), ncol=ncol(testCounts))

#make test spUnsupervised and spSwarm
cObjSng <- spCounts(testCounts[,s2], testErcc[,s2])
cObjMul <- spCounts(testCounts[,!s2], testErcc[,!s2])
testUns <- spUnsupervised(cObjSng, max=250, max_iter=1000)

testSwa <- spSwarm(
    cObjMul,
    testUns,
    distFun=distToSlice,
    maxiter=100,
    swarmsize=500,
    cores=2
)

#save
save(testErcc, file="data/testErcc.rda", compress="bzip2")
save(testCounts, file="data/testCounts.rda", compress="bzip2")
save(testUns, file="data/testUns.rda", compress="bzip2")
save(testSwa, file="data/testSwa.rda", compress="bzip2")