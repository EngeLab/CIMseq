data(expData)
data(unsupervised)

classification <- getData(unsupervised, "mclust")$classification
counts.log <- getData(expData, "counts.log")
sampleType <- getData(expData, "sampleType")
dbl <- counts.log[ ,sampleType == "Doublet"]
sng <- counts.log[ ,sampleType == "Singlet"]



.averageGroupExpression <- function(classes, sng) {
    classes <- unique(classes)
    means <- lapply(classes, function(x) {
        ingroup <- classes == x
        log2(rowMeans(2^sng[,ingroup]))
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- classes
    return(means)
}

clusterMeans <- .averageGroupExpression(classification, sng)
fractions <- rep(1.0/(dim(clusterMeans)[2]), (dim(clusterMeans)[2]))





maxs <- order(apply(counts.log, 1, max), decreasing=T)
clusterMeans.top2000 <- 2^clusterMeans[maxs[1:2000],]
dbl.top2000 <- 2^dbl[maxs[1:2000],]


optim.res.all.2 <- lapply(1:(dim(dbl)[2]), function(i) {
    optimx(par=fractions, fn=dist.to.slice, gr=NULL, cell.types=clusterMeans.top2000, slice=dbl.top2000[,i], method=c("L-BFGS-B"), lower=0.0, upper=1.0)
})


groups.mat <- as.matrix(data.frame(lapply(optim.res.all.2, function(x) {unlist(x[1:12])})))
heatmap(t(groups.mat), col=colorRampPalette(c("white", "red"))(100), ColSideColors=jet.colors(max(mod1$classification))[unique(mod1$classification)])

dist.to.slice <- function(fractions, cell.types, slice) {
    a <- make.synthetic.slice(cell.types, fractions)
    if(any(is.na(a))) {
        cat("NA in make synthetic slice!\n")
        cat(paste(fractions, sep="\t"), "\n")
    }
    cost <- sum(abs(a - slice))
    cost
}

# Distance function for the optimization - could possibly be better.
dist.to.slice <- function(fractions, cell.types, slice) {
    
    tmp <- make.synthetic.slice(cell.types, fractions)
    fractions <- tmp[[1]]
    a <- tmp[[2]]
    if(any(is.na(a))) {
        cat("NA in make synthetic slice!\n")
        cat(paste(fractions, sep="\t"), "\n")
    }
    cost <- sum(abs(a - slice))
    cost
}

# Support function for the optimization
make.synthetic.slice <- function(cell.types, fractions) {
    fractions <- fractions/sum(fractions)
    res <- apply(cell.types, 1, function(x) {sum(x*fractions)})
    return(list(fractions, res))
}

##############try with test data
classification <- getData(unsupervised, "mclust")$classification
counts.log <- getData(expData, "counts.log")
sampleType <- getData(expData, "sampleType")
sng <- counts.log[ ,sampleType == "Singlet"]

clusterMeans <- .averageGroupExpression(classification, sng)
clusterMeans <- clusterMeans[ ,c('4', '7', '8')]
fractions <- rep(1.0/(dim(clusterMeans)[2]), (dim(clusterMeans)[2]))

maxs <- order(apply(counts.log, 1, max), decreasing=T)
clusterMeans.top2000 <- 2^clusterMeans[maxs[1:2000],]
dbl.top2000 <- 2^testData[maxs[1:2000],]


optim.res.all.2 <- lapply(1:(dim(testData)[2]), function(i) {
    optimx(par=fractions, fn=dist.to.slice, gr=NULL, cell.types=clusterMeans.top2000, slice=dbl.top2000[,i], method=c("L-BFGS-B"), lower=0, upper=1.0)
})


groups.mat <- as.matrix(data.frame(lapply(optim.res.all.2, function(x) {unlist(x[1:3])})))
colnames(groups.mat) <- colnames(testData)









##################try slsqp


library(nloptr)
#, lower=0, upper=1.0
slsqp <- lapply(1:(dim(testData)[2]), function(i) {
    slsqp(x0=fractions, fn=dist.to.slice, gr=NULL, cell.types=clusterMeans.top2000, slice=dbl.top2000[,i])
})







