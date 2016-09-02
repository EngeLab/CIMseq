
.averageGroupExpression <- function(classes, sng) {
    c <- unique(classes)
    means <- lapply(c, function(x) {
        rowMeans(sng[,classes == x])
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
}


assembleTestData <- function() {
    load('data/expData.rda')
    load('data/unsupervised.rda')
    
    counts <- getData(expData, "counts")
    sampleType <- getData(expData, "sampleType")
    
    maxs <- order(apply(counts.log, 1, max), decreasing=T)
    counts <- counts[maxs[1:2000], ]
    sng <- counts[ ,sampleType == "Singlet"]
    
    classification <- getData(unsupervised, "mclust")$classification
    means <- .averageGroupExpression(classification, sng)
    
    #doublets
    comb <- combn(unique(classification), 2)
    
    for( i in 1:ncol(comb)) {
        currentMult <- rowMeans(data.frame(means[ ,comb[1,i]], means[ ,comb[2,i]]))
        name <- paste(comb[1,i], comb[2,i], sep=".")
        
        if( i == 1) {
            dataset <- data.frame(one = currentMult)
            names <- name
        } else {
            dataset <- cbind(dataset, currentMult)
            names <- c(names, name)
        }
    }
    
    #triplets
    comb <- combn(unique(classification), 3)
    
    for( i in 1:ncol(comb)) {
        currentMult <- rowMeans(data.frame(means[ ,comb[1,i]], means[ ,comb[2,i]], means[ ,comb[3,i]]))
        name <- paste(comb[1,i], comb[2,i], comb[3,i], sep=".")
        
        dataset <- cbind(dataset, currentMult)
        names <- c(names, name)
    }
    
    colnames(dataset) <- names
    testData <- dataset
    save(testData, file="data/testData.rda", compress="bzip2")
}


syntheticTestData <- function(save=FALSE) {
    singlets <- .syntheticSinglets()
    multuplets <- .syntheticMultuplets()
    names <- c(paste("s", colnames(singlets), sep="."), paste("m", colnames(multuplets), sep="."))
    
    syntheticData <- cbind(singlets, multuplets)
    colnames(syntheticData) <- names
    syntheticData <- as.matrix(syntheticData)
    
    if(save == TRUE) {
        save(syntheticData, file="data/syntheticData.rda", compress="bzip2")
    }
    
    return(syntheticData)
}

.syntheticSinglets <- function(save = FALSE) {
    ngenes <- 2000
    ncells <- 100
    cellTypes <- 10
    for( i in 1:cellTypes) {
        set.seed(i)
        meanExprs <- 2^runif(ngenes, 0, 5)
        counts <- matrix(rnbinom(ngenes*ncells, mu=meanExprs, size=i), nrow=ngenes)
        if( i == 1 ) {
            singlets <- counts
        } else {
            singlets <- cbind(singlets, counts)
        }
    }
    colnames(singlets) <- paste(sort(rep(letters, ncells))[1:(cellTypes*ncells)], 1:100, sep="")
    singlets <- as.data.frame(singlets)
    
    return(singlets)
}

.syntheticMultuplets <- function(save = FALSE) {
    singlets <- .syntheticSinglets()
    cObj <- spCounts(as.matrix(singlets), counts.ercc=matrix(), sampleType="[A-Z]")
    uObj <- spUnsupervised(cObj, max=1000, max_iter = 10)
    colnames(singlets) <- getData(uObj, "mclust")$classification
    
    mean <- getData(uObj, "groupMeans")
    
    #doublets
    combos <- combn(unique(colnames(singlets)), 2)
    
    for(y in 1:ncol(combos)) {
        current <- combos[ ,y]
        new <- data.frame(rowMeans(mean[ , colnames(mean) %in% current]))
        
        if( y == 1 ) {
            multuplets <- new
            names <- paste(current, collapse="")
        } else {
            multuplets <- cbind(multuplets, new)
            names <- c(names, paste(current, collapse=""))
        }
    }
    
    colnames(multuplets) <- names
    
    #triplets
    newNames <- colnames(multuplets)
    combos <- combn(unique(colnames(singlets)), 3)
    
    for(u in 1:ncol(combos)) {
        current <- combos[ ,u]
        new <- data.frame(rowMeans(mean[ , colnames(mean) %in% current]))
        multuplets <- cbind(multuplets, new)
        newNames <- c(newNames,  paste(current, collapse=""))
    }
    colnames(multuplets) <- newNames
    
    return(multuplets)
}















