
#'@include All-classes.R
NULL

#' spSwarm
#'
#' Subtitle
#'
#' Description
#'
#' @name spSwarm
#' @rdname spSwarm
#' @aliases spSwarm
#' @param spUnsupervised An spUnsupervised object.
#' @param limit Randomly select a limited number of samples to perform swarm optimization on.
#' @param maxiter pySwarm argument indicating the maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of particals in the swarm.
#' @param minstep pySwarm argument indicating the stepsize of swarm’s best position before search termination.
#' @param minfunc pySwarm argument indicating the minimum change of swarm’s best objective value before search termination.
#' @param cutoff The cutoff used to generate codedSwarm.
#' @param cores The number of cores to be used while running spSwarm.
#' @param n Data to extract from spSwarm object.
#' @param .Object Internal object.
#' @param object spSwarm object.
#' @param arguments Argumetns passed to the spSwarm function.
#' @param spSwarm The spSwarm results.
#' @param codedSwarm The spSwarm results after the cutoff has been applied.
#' @param x A spSwarm object.
#' @param ... additional arguments to pass on
#' @return spSwarm output.
#' @author Jason T. Serviss
#' @keywords spSwarm, pySwarm
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname spSwarm
#' @export

setGeneric("spSwarm", function(spUnsupervised, ...
){ standardGeneric("spSwarm") })

#' @importFrom rPython python.exec python.assign python.get
#' @importFrom doMC registerDoMC
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @rdname spSwarm
#' @export
setMethod("spSwarm", "spUnsupervised",
function(
    spUnsupervised,
    limit = "none",
    maxiter = 10,
    swarmsize = 150,
    minstep = 1e-16,
    minfunc = 1e-16,
    cutoff = 0.2,
    cores=1,
    ...
){
    
    #input and input checks
    counts <- getData(spUnsupervised, "counts")
    sampleType <- getData(spUnsupervised, "sampleType")
    selectInd <- getData(spUnsupervised, "selectInd")
    
    #calculate average expression
    groupMeans <- getData(spUnsupervised, "groupMeans")
    
    #calculate fractions
    fractions <- rep(1.0/(dim(groupMeans)[2]), (dim(groupMeans)[2]))
    
    #subset top genes for use with optimization
    cellTypes <- as.data.frame(groupMeans[selectInd, ])
    slice <- as.data.frame(counts[selectInd, sampleType == "Multuplet"]) #note there is an error here if there is only one multuplet due to the fact that the subset resluts in a data.frame with a colname that does not match the colnames in the counts variable
    
    #add index for reordering in python
    cellTypes$index <- 1:nrow(cellTypes)
    slice$index <- 1:nrow(slice)
    
    ##run pySwarm
    .defineImport(cellTypes, slice, fractions)
    .defineCost()
    .con1()
    .definePySwarm()
    result <- .runPySwarm(
        cellTypes,
        slice,
        fractions,
        limit,
        maxiter,
        swarmsize,
        minstep,
        minfunc,
        cores
    )
        
    #process and return results
    finalResult <- .processResults(result)
    encodedResult <- .multiHOTencoding(finalResult, spCounts, cutoff)
    
    #create object
    new("spSwarm",
        spSwarm=finalResult,
        codedSwarm=encodedResult,
        spUnsupervised=spUnsupervised,
        arguments = list(
            maxiter=maxiter,
            swarmsize=swarmsize,
            minstep=minstep,
            minfunc=minfunc
        )
    )
})

#############################################
#                                           #
#              Unit Functions               #
#                                           #
#############################################

##define optimization and constraint functions
##import data to python session
.defineImport <- function(cellTypes, slice, fractions) {
    cmd1 <- 'import pandas as pd'
    cmd2 <- 'import numpy as np'
    cmd3 <- 'import math'
    python.exec(cmd1)
    python.exec(cmd2)
    python.exec(cmd3)
    
    python.assign('cellTypesDictionary', cellTypes)
    python.assign('sliceDictionary', slice)
    python.assign('fractions', fractions)
    
    cmd4 <- 'cellTypes = pd.DataFrame(cellTypesDictionary)'
    cmd5 <- 'slice = pd.DataFrame(sliceDictionary)'
    cmd6 <- 'fractions = np.asarray(fractions)'
    python.exec(cmd4)
    python.exec(cmd5)
    python.exec(cmd6)
    
    cmd7 <- 'cellTypes.sort_values(by=\'index\', ascending=\'True\')'
    cmd8 <- 'slice.sort_values(by=\'index\', ascending=\'True\')'
    python.exec(cmd7)
    python.exec(cmd8)
    
    cmd9 <- 'cellTypes = cellTypes.drop(\'index\', 1)'
    cmd10 <- 'slice = slice.drop(\'index\', 1)'
    python.exec(cmd9)
    python.exec(cmd10)
    
    return(list(
        cmd1,
        cmd2,
        cmd3,
        cmd4,
        cmd5,
        cmd6,
        cmd7,
        cmd8,
        cmd9,
        cmd10
    ))
}

##define cost function(s) in python session
.defineCost <- function() {
    
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
    return(list(cmd1, cmd2))
}

#constraint 1
.con1 <- function() {
    cmd1 <- 'def con1(fractions, *args):
        if sum(fractions) > 0.1:
            return 0
        else:
            return -1'
    
    python.exec(cmd1)
    return(cmd1)
}

#optimization
.definePySwarm <- function() {
    python.exec('from pyswarm import pso')
    
    cmd1 <- 'def optimize(cellTypes, slice, fractions, col, maxiter, swarmsize, minstep, minfunc):
        lb = np.asarray([0] * len(fractions))
        ub = np.asarray([1] * len(fractions))
        name = slice.columns.values[col]
        args = (cellTypes, slice, col)
        xopt, fopt = pso(
            distToSlice,
            lb,
            ub,
            args=args,
            f_ieqcons=con1,
            maxiter=maxiter,
            swarmsize=swarmsize,
            minstep=minstep,
            minfunc=minfunc
        )
        dictionary = dict(zip(list(cellTypes.columns.values), xopt.tolist()))
        return { \'xopt\':dictionary, \'fopt\':fopt, \'name\':name }'
        
    python.exec(cmd1)
    return(cmd1)
}

##run optimization
.runPySwarm <- function(
    cellTypes,
    slice,
    fractions,
    limit,
    maxiter,
    swarmsize,
    minstep,
    minfunc,
    cores
){
    result <- list()
    
    ##setup parallel processing
    registerDoMC(cores)
    
    if(limit == "none") {
        top <- 0:(ncol(slice) - 2)
    } else {
        top <- sample(
            0:(ncol(slice) - 2),
            limit,
            replace=FALSE
        )
    }
    
    
    while(TRUE) {
        
        print(paste(length(top), " multuplets left to analyze.", sep=""))

        loopOutput <- foreach(i = 1:cores, top=top, .combine = append) %dopar% {
            python.exec(
                paste(
                    'result = optimize(cellTypes, slice, fractions, col=',
                    top,
                    ', maxiter=', maxiter,
                    ', swarmsize=', swarmsize,
                    ', minstep=', minstep,
                    ', minfunc=', minfunc,
                    ')',
                    sep=""
                )
            )
            result <- list(
                currentFopt = python.get(paste('result[\'fopt\']')),
                currentXopt = python.get(paste('result[\'xopt\']')),
                name = python.get(paste('result[\'name\']'))
            )
            
        }
        
        result <- append(result, loopOutput)
        top <- top[(1:cores)*-1]
        if(length(top) == 0) {break}
        
    }
    
    registerDoSEQ()
    return(result)
}



.processResults <- function(result) {
    
    #extract xopt
    xopt <- data.frame(t(data.frame(result[which(names(result) == "currentXopt")])))
    xopt <- xopt[, order(colnames(xopt))]
    rownames(xopt) <- 1:nrow(xopt)
    
    #normalize xopt with sum
    xopt <- xopt/rowSums(xopt)
    
    #extract cost (fopt)
    xopt$fopt <- unlist(result[which(names(result) == "currentFopt")])
    
    #extract names (fopt)
    xopt$sampleName <- unlist(result[which(names(result) == "name")])
    
    #add sample and fopt data
    names <- c("sampleName", "fopt")
    finalRes <- cbind(subset(xopt, select=names), xopt[ ,!colnames(xopt) %in% names])
    return(finalRes)
}

.multiHOTencoding <- function(optResult, spCounts, cutoff) {
    
    uu <- c("names", "fopt")
    hold <- data.frame(fopt = optResult[ , colnames(optResult) %in% uu])
    x <- optResult[ , !colnames(optResult) %in% uu]
    
    if(class(cutoff) == "numeric") {
        
        for(p in 1:nrow(x)) {
            x[p,][x[p,] < cutoff] <- 0
        }
        
    } else {
        counts <- getData(spCounts, "counts")
        counts.ercc <- getData(spCounts, "counts.ercc")
        sampleType <- getData(spCounts, "sampleType")
        
        frac.ercc <- 100 * (colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts)))
        frac.ercc <- frac.ercc[sampleType == "Multuplet"]
        
        for(p in 1:nrow(x)) {
            cutoff <- frac.ercc[p] / ncol(x)
            x[p,][x[p,] < cutoff] <- 0
        }
        
    }
    
    return(cbind(hold, x))
}

#.PCAsubset <- function(counts, sampleType, cutoff) {

    #subset singlets
    #    sng <- counts[ ,sampleType == "Singlet"]
    
    #remove 0's
    #    c <- counts[rowMeans(sng) != 0, ]
    
    #run PCA
    #    pca <- prcomp(t(c), center=TRUE, scale=TRUE)
    
    #order and return genes
    #    rotation <- as.data.frame(pca$rotation)
    #    genes <- rownames(rotation[order(rotation$PC1), ])
    #    cut <- cutoff/2
    #    selected <- c(genes[1:cut], genes[(length(genes) - (cut-1)):length(genes)])
    #    return(selected)
    #}

#' changeCutoff
#'
#' Subtitle
#'
#' Description
#'
#' @name changeCutoff
#' @rdname changeCutoff
#' @aliases changeCutoff
#' @param spSwarm An spSwarm object.
#' @param cutoff The fraction below which a connection should not be considered.
#' @param ... additional arguments to pass on
#' @return spSwarm output.
#' @author Jason T. Serviss
#' @keywords changeCutoff
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname changeCutoff
#' @export

setGeneric("changeCutoff", function(spSwarm, ...
){ standardGeneric("changeCutoff") })

#' @rdname changeCutoff
#' @export

setMethod("changeCutoff", "spSwarm",
function(
    spSwarm,
    cutoff,
    ...
){
    spCounts <- getData(spSwarm, "spCounts")
    swarmResults <- getData(spSwarm, "spSwarm")
    
    uu <- c("names", "fopt")
    
    hold <- swarmResults[ , colnames(swarmResults) %in% uu]
    x <- swarmResults[ , !colnames(swarmResults) %in% uu]
    
    if(class(cutoff) == "numeric") {
        
        for(p in 1:nrow(x)) {
            x[p,][x[p,] < cutoff] <- 0
        }
        
    } else {
        counts <- getData(spCounts, "counts")
        counts.ercc <- getData(spCounts, "counts.ercc")
        sampleType <- getData(spCounts, "sampleType")
        
        frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
        medianErccSng <- median(frac.ercc[sampleType == "Singlet"])
        extimatedCells <- medianErccSng/frac.ercc[sampleType == "Multuplet"]
        
        for(p in 1:nrow(x)) {
            cutoff <- extimatedCells[p] / ncol(x)
            print(cutoff)
            x[p,][x[p,] < cutoff] <- 0
        }
        
    }
    
    codedSwarm <- cbind(hold, x)
    spSwarm@codedSwarm <- codedSwarm
    return(spSwarm)
})


