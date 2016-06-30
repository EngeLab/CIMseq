
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
    ...
){
    #input and input checks
    classification <- getData(spUnsupervised, "mclust")$classification
    counts.log <- getData(spCounts, "counts.log")
    sampleType <- getData(spCounts, "sampleType")
    sng <- counts.log[ ,sampleType == "Singlet"]
    
    #calculate average expression
    clusterMeans <- .averageGroupExpression(classification, sng)
    
    #calculate fractions
    fractions <- rep(1.0/(dim(clusterMeans)[2]), (dim(clusterMeans)[2]))
    
    #subset top 2000 genes for use with optimization
    #maxs <- order(apply(counts.log, 1, max), decreasing=T)
    maxs <- order(apply(counts.log, 1, var), decreasing=T)
    cellTypes <- 2^clusterMeans[maxs[1:2000],]
    slice <- 2^testData2[maxs[1:2000],]
    
    ##run pySwarm
    .defineImport(cellTypes, slice, fractions)
    .defineCost()
    .con1()
    .con2()
    .constraints()
    .definePySwarm()
    result <- .runPySwarm(cellTypes, slice, fractions)
    
    return(result)
})

##define optimization and constraint functions
#constraint 1
.con1 <- function() {
    cmd1 <- 'def con1(fractions, *args):
        if sum(fractions) == 1:
            return 0
        else:
            return -1'
    
    python.exec(cmd1)
}

#constraint 2
.con2 <- function() {
    cmd1 <- 'from collections import Counter'
    cmd2 <- 'def con2(fractions, *args):
        limit = 1/np.float64(len(fractions))
        oc1 = list(set([i for i in fractions if i < limit and i != 0.0]))
        if sum([Counter(fractions)[k] for k in oc1]) == 0:
            return 0
        else:
            return -1'
    
    python.exec(cmd1)
    python.exec(cmd2)
}

#all constraints
.constraints <- function() {
    cmd1 <- 'def constraints(fractions, *args):
        cons1 = con1(fractions, *args)
        cons2 = con2(fractions, *args)
        return [cons1, cons2]'
    
    python.exec(cmd1)
}

#optimization
.definePySwarm <- function() {
    python.exec('from pyswarm import pso')
    
    cmd1 <- 'def optimize(cellTypes, slice, fractions, col):
        lb = np.asarray([0] * len(fractions))
        ub = np.asarray([1] * len(fractions))
        args = (cellTypes, slice, col)
        xopt, fopt = pso(distToSlice, lb, ub, args=args, f_ieqcons=constraints)
        dictionary = dict(zip(list(cellTypes.columns.values), xopt.tolist))
        return { \'xopt\':dictionary, \'fopt\':fopt }'
        
    #python.exec(cmd1)
    python.exec(cmd1)
}

##run optimization
.runPySwarm <- function(cellTypes, slice, fractions) {
    result <- list()
    for(pp in 0:(ncol(slice) - 1)) {
        print(paste("analyzing multuplet number: ", pp, sep=""))
        python.exec(paste('result = optimize(cellTypes, slice, fractions, col=', pp, ')', sep=""))
        result[[(pp+1)]] = list(
            currentFopt = python.get(paste('result[\'fopt\']')),
            currentXopt = python.get(paste('result[\'xopt\']'))
        )
    }
    return(result)
}


#############################################
#                                           #
#           General Functions               #
#                                           #
#############################################

##calculates the average expression for each singlet group
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

##import data to python session
.defineImport <- function(cellTypes, slice, fractions) {
    python.exec('import pandas as pd')
    python.exec('import numpy as np')
    
    python.assign('cellTypesDictionary', as.data.frame(cellTypes))
    python.assign('sliceDictionary', as.data.frame(slice))
    python.assign('fractions', fractions)
    
    python.exec('cellTypes = pd.DataFrame(cellTypesDictionary)')
    python.exec('slice = pd.DataFrame(sliceDictionary)')
    python.exec('fractions = np.asarray(fractions)')
}

##define cost function(s) in python session
.defineCost <- function() {
    python.exec('import pandas as pd')
    python.exec('import numpy as np')
    python.exec('import math')
    
    cmd1 <- 'def makeSyntheticSlice(cellTypes, fractions):
    #s = sum(fractions)
    
    #for i, f in enumerate(fractions):
    #fractions[i] = np.float64(fractions[i]) / s
    
        func = lambda x: sum(np.asarray(x) * np.asarray(fractions))
        return cellTypes.apply(func, axis=1)'
    
    cmd2 <- 'def distToSlice(fractions, *args):
        cellTypes, slice, col = args
        a = makeSyntheticSlice(cellTypes, fractions)
        for i in a:
            if math.isnan(i):
                print "NaN in make synthetic slice!"
                print a
    
        diff = []
        for index, row in cellTypes.iterrows():
            d = a[index] - slice.iloc[index, col]
            diff.append(abs(d))
        cost = sum(diff)
        return cost'
    
    python.exec(cmd1)
    python.exec(cmd2)
}






