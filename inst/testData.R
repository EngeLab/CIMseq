
assembleTestData <- function() {
    load('data/expData.rda')
    load('data/unsupervised.rda')
    
    counts.log <- getData(expData, "counts.log")
    sampleType <- getData(expData, "sampleType")
    sng <- counts.log[ ,sampleType == "Singlet"]
    classification <- getData(unsupervised, "mclust")$classification
    
    testData <- data.frame()

    ##4 and 7 as doublets
    g1 <- sng[ ,classification == 4]
    g2 <- sng[ ,classification == 7]
    int <- min(ncol(g1), ncol(g2))
    
    for(oo in 1:int) {
        a <- g1[ ,oo]
        for(tt in 1:int) {
            b <- g2[ ,tt]
            m <- data.frame(rowMeans(cbind(a,b)))
            
            if(tt == 1 & oo == 1) {
                testData <- m
            } else {
                testData <- cbind(testData, m)
            }
        }
    }
    names <- paste("fourSeven", 1:(int^2), sep="")
    
    ##4 and 9 as doublets
    g1 <- sng[ ,classification == 4]
    g2 <- sng[ ,classification == 9]
    int <- min(ncol(g1), ncol(g2))
    
    for(oo in 1:int) {
        a <- g1[ ,oo]
        for(tt in 1:int) {
            b <- g2[ ,tt]
            m <- data.frame(rowMeans(cbind(a,b)))
            testData <- cbind(testData, m)
        }
    }
    names <- c(names, paste("fourNine", 1:(int^2), sep=""))
    
    ##7 and 9 as doublets
    g1 <- sng[ ,classification == 7]
    g2 <- sng[ ,classification == 9]
    int <- min(ncol(g1), ncol(g2))
    
    for(oo in 1:int) {
        a <- g1[ ,oo]
        for(tt in 1:int) {
            b <- g2[ ,tt]
            m <- data.frame(rowMeans(cbind(a,b)))
            testData <- cbind(testData, m)
        }
    }
    names <- c(names, paste("sevenNine", 1:(int^2), sep=""))
    
    ##4, 7, and 9 as triplets
    g1 <- sng[ ,classification == 4]
    g2 <- sng[ ,classification == 7]
    g3 <- sng[ ,classification == 9]

    int <- min(ncol(g1), ncol(g2), ncol(g3))
    
    for(oo in 1:int) {
        a <- g1[ ,oo]
        for(tt in 1:int) {
            b <- g2[ ,tt]
            for(rr in 1:int) {
                c <- g3[ ,rr]
                m <- data.frame(rowMeans(cbind(a,b,c)))
                testData <- cbind(testData, m)
            }
        }
    }
    names <- c(names, paste("fourSevenNine", 1:(int^3), sep=""))
    colnames(testData) <- names
    
    save(testData, file="data/testData.rda", compress="bzip2")
    return(testData)
}




.averageGroupExpression <- function(classes, sng) {
    c <- unique(classes)
    means <- lapply(c, function(x) {
        rowMeans(sng[,classes == x])
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
}


assembleTestData2 <- function() {
    load('data/expData.rda')
    load('data/unsupervised.rda')
    
    counts <- getData(expData, "counts")
    sampleType <- getData(expData, "sampleType")
    sng <- counts[ ,sampleType == "Singlet"]
    classification <- getData(unsupervised, "mclust")$classification
    means <- .averageGroupExpression(classification, sng)
    
    #doublets
    comb <- combn(seq(1, length(unique(classification)), 1), 2)
    
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
    comb <- combn(seq(1, length(unique(classification)), 1), 3)
    
    for( i in 1:ncol(comb)) {
        currentMult <- rowMeans(data.frame(means[ ,comb[1,i]], means[ ,comb[2,i]], means[ ,comb[3,i]]))
        name <- paste(comb[1,i], comb[2,i], comb[3,i], sep=".")
        
        dataset <- cbind(dataset, currentMult)
        names <- c(names, name)
    }
    
    colnames(dataset) <- names
    testData2 <- dataset
    save(testData2, file="data/testData2.rda", compress="bzip2")
}


syntheticSinglets <- function(save = FALSE) {
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
    
    if( save == TRUE ) {
        save(singlets, file="data/syntheticSng.rda", compress="bzip2")
    }
    
    return(singlets)
}

syntheticMultupletsA <- function(save = FALSE) {
    singlets <- syntheticSinglets()
    newColNames <- unlist(strsplit(colnames(singlets), "[0-9]"))
    colnames(singlets) <- newColNames[newColNames != ""]
    
    mean <- .averageGroupExpression(classes = colnames(singlets), sng = singlets)
    
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
    
    if( save == TRUE ) {
        save(multuplets, file="data/syntheticMultupletsA.rda", compress="bzip2")
    }
    
    return(multuplets)
}

#syntheticMultupletsB <- function() {
#    singlets <- syntheticSinglets()
#    names <- unique(colnames(singlets))
#
#    for( i in 1:length(names) ) {
#        type1 <- singlets[ ,colnames(singlets) == names[i]]
#
#        for(y in 1:length(names) ) {
#            type2 <- singlets[ ,colnames(singlets) == names[y]]
#
#            for( j in 1:ncol(type1) ) {
#
#                for( k in 1:ncol(type2) ) {
#                    cell1 <- type1[ ,j]
#                    cell2 <- type2[ ,k]
#                    name1 <- paste(names[i], j, sep="")
#                    name2 <- paste(names[y], k, sep="")
#                    comb <- data.frame(
#                        rowMeans(
#                            cbind(
#                                data.frame(cell1),
#                                data.frame(cell2)
#                            )
#                        )
#                    )
#
#                    if( i ==1 & y==1 & j==1 & k==1) {
#                        output <- comb
#                        newNames <- paste(name1, name2, sep="")
#                    } else {
#                        output <- cbind(output, comb)
#                        newNames <- c(newNames, paste(name1, name2, sep=""))
#                    }
#
#                }
#            }
#        }
#    }
#}



















