
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
    
    fourSeven <- rowMeans(data.frame(means[ ,'4'], means[ ,'7']))
    sevenNine <- rowMeans(data.frame(means[ ,'7'], means[ ,'9']))
    fourNine <- rowMeans(data.frame(means[ ,'4'], means[ ,'9']))
    fourSevenNine <- rowMeans(data.frame(means[ ,'4'], means[ ,'7'], means[ ,'9']))
    
    testData2 <- data.frame(fourSeven = fourSeven, sevenNine = sevenNine, fourNine = fourNine, fourSevenNine = fourSevenNine)
    save(testData2, file="data/testData2.rda", compress="bzip2")
}























