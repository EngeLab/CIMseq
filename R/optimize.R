
#'@include All-classes.R
NULL

#' pySwarm
#'
#' Subtitle
#'
#' Description
#'
#' @name pySwarm
#' @rdname pySwarm
#' @aliases pySwarm
#' @param spCounts An spCounts object.
#' @param spUnsupervised An spUnsupervised object.
#' @param ... additional arguments to pass on
#' @return PySwarm output.
#' @author Jason T. Serviss
#' @keywords pySwarm
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname pySwarm
#' @export

setGeneric("pySwarm", function(spCounts, ...
){ standardGeneric("pySwarm") })

#' @importFrom rPython python.exec python.assign python.get
#' @rdname pySwarm
#' @export
setMethod("pySwarm", "spCounts",
function(
    spCounts,
    spUnsupervised,
    limit = "none",
    ...
){
    
    #input and input checks
    counts <- getData(spCounts, "counts")
    sampleType <- getData(spCounts, "sampleType")
    sng <- counts[ ,sampleType == "Singlet"]
    mult <- counts[ ,sampleType == "Multuplet"]
    classification <- getData(spUnsupervised, "mclust")$classification
    
    #calculate average expression
    clusterMeans <- .averageGroupExpression(classification, sng)
    
    #calculate fractions
    fractions <- rep(1.0/(dim(clusterMeans)[2]), (dim(clusterMeans)[2]))
    
    #subset top 2000 genes for use with optimization
        #use max expression
        maxs <- order(apply(counts, 1, max), decreasing=T)
        cellTypes <- as.data.frame(clusterMeans[maxs[1:2000],]) #log scale
        slice <- as.data.frame(mult[maxs[1:2000],]) #log scale
        
        #use variance
        #maxs <- order(apply(counts, 1, var), decreasing=T)
        #cellTypes <- as.data.frame(clusterMeans[maxs[1:2000],]) #log scale
        #slice <- as.data.frame(testData2[maxs[1:2000],]) #log scale
        
        #use PCA
        #genes <- .PCAsubset(counts, sampleType, 100)
        #cellTypes <- as.data.frame(clusterMeans[rownames(clusterMeans) %in% genes, ])
        #slice <- as.data.frame(testData2[rownames(clusterMeans) %in% genes, ])
    
    #add index for reordering in python
    cellTypes$index <- 1:nrow(cellTypes)
    slice$index <- 1:nrow(slice)
    
    ##run pySwarm
    .defineImport(cellTypes, slice, fractions)
    .defineCost()
    .con1()
    .definePySwarm()
    result <- .runPySwarm(cellTypes, slice, fractions, limit)
    
    return(result)
})

##define optimization and constraint functions
#constraint 1
.con1 <- function() {
    cmd1 <- 'def con1(fractions, *args):
        if sum(fractions) > 0.1:
            return 0
        else:
            return -1'
    
    python.exec(cmd1)
}

#optimization
.definePySwarm <- function() {
    python.exec('from pyswarm import pso')
    
    cmd1 <- 'def optimize(cellTypes, slice, fractions, col):
        lb = np.asarray([0] * len(fractions))
        ub = np.asarray([1] * len(fractions))
        args = (cellTypes, slice, col)
        xopt, fopt = pso(distToSlice, lb, ub, args=args, f_ieqcons=con1, maxiter=10, swarmsize=250, minstep=1e-16, minfunc=1e-16)
        dictionary = dict(zip(list(cellTypes.columns.values), xopt.tolist()))
        return { \'xopt\':dictionary, \'fopt\':fopt }'
        
    python.exec(cmd1)
}

##run optimization
.runPySwarm <- function(cellTypes, slice, fractions, limit) {
    result <- list()
    
    if( limit == "none") {
        top <- 0:(ncol(slice) - 1)
    } else {
        top <- sample(
            0:(ncol(slice) - 1),
            limit,
            replace=FALSE
        )
    }
    
    for(pp in 1:length(top)) {
        print(paste("analyzing multuplet number: ", pp, sep=""))
        col <- top[pp]
        python.exec(paste('result = optimize(cellTypes, slice, fractions, col=', col, ')', sep=""))
        result[[(pp)]] = list(
            currentFopt = python.get(paste('result[\'fopt\']')),
            currentXopt = python.get(paste('result[\'xopt\']'))
        )
    }
    names(result) <- colnames(slice)[top]
    return(result)
}


#############################################
#                                           #
#           General Functions               #
#                                           #
#############################################

##calculates the average expression for each singlet group
.averageGroupExpression <- function(classes, sng) {
    c <- unique(classes)
    means <- lapply(c, function(x) {
        rowMeans(sng[,classes == x])
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
}

##import data to python session
.defineImport <- function(cellTypes, slice, fractions) {
    python.exec('import pandas as pd')
    python.exec('import numpy as np')
    
    python.assign('cellTypesDictionary', cellTypes)
    python.assign('sliceDictionary', slice)
    python.assign('fractions', fractions)
    
    python.exec('cellTypes = pd.DataFrame(cellTypesDictionary)')
    python.exec('slice = pd.DataFrame(sliceDictionary)')
    python.exec('fractions = np.asarray(fractions)')
    
    python.exec('cellTypes.sort_values(by=\'index\', ascending=\'True\')')
    python.exec('slice.sort_values(by=\'index\', ascending=\'True\')')

    python.exec('cellTypes = cellTypes.drop(\'index\', 1)')
    python.exec('slice = slice.drop(\'index\', 1)')
}

##define cost function(s) in python session
.defineCost <- function() {
    python.exec('import pandas as pd')
    python.exec('import numpy as np')
    python.exec('import math')
    
    cmd1 <- 'def makeSyntheticSlice(cellTypes, fractions):
        func = lambda x: sum(np.asarray(x) * np.asarray(fractions))
        return cellTypes.apply(func, axis=1)'
    
    cmd2 <- 'def distToSlice(fractions, *args):
    cellTypes, slice, col = args
    normFractions = fractions / sum(fractions)
    a = makeSyntheticSlice(cellTypes, normFractions)
    for i in a:
        if math.isnan(i):
            print "NaN in make synthetic slice!"
    
    diff = []
    for index in range(cellTypes.shape[0]):
        d = a.iloc[index] - slice.iloc[index, col]
        diff.append(abs(d))
    cost = sum(diff)
    return cost'
    
    python.exec(cmd1)
    python.exec(cmd2)
}

.PCAsubset <- function(counts, sampleType, cutoff) {
    
    #subset singlets
    sng <- counts[ ,sampleType == "Singlet"]
    
    #remove 0's
    c <- counts[rowMeans(sng) != 0, ]
    
    #run PCA
    pca <- prcomp(t(c), center=TRUE, scale=TRUE)
    
    #order and return genes
    rotation <- as.data.frame(pca$rotation)
    genes <- rownames(rotation[order(rotation$PC1), ])
    cut <- cutoff/2
    selected <- c(genes[1:cut], genes[(length(genes) - (cut-1)):length(genes)])
    return(selected)
}




