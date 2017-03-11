#'@include All-classes.R
NULL

# Helper functions to classify count data into clusters

#' spSwarm
#'
#' Subtitle
#'
#' Description
#'
#' @name spSwarm
#' @rdname spSwarm
#' @aliases spSwarm
#' @param spCounts an spCount object.
#' @param maxiter pySwarm argument indicating the maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of particals in the swarm.
#' @param cores The number of cores to be used while running spRSwarm.
#' @param arguments Argumetns passed to the spRSwarm function.
#' @param spRSwarm The spRSwarm results.
#' @param object spRSwarm object.
#' @param n Data to extract from spRSwarm object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spRSwarm output.
#' @author Jason T. Serviss
#' @keywords spSwarm
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

setGeneric("spSwarm", function(spCounts, spUnsupervised, ...
){ standardGeneric("spSwarm") })


#' @importFrom parallel mclapply
#' @importFrom pso psoptim
#' @rdname spSwarm
#' @export

setMethod("spSwarm", c("spCounts", "spUnsupervised"), function(
    spCounts,
    spUnsupervised,
    distFun = distToSlice,
    selectInd = TRUE,
    maxiter = 10,
    swarmsize = 150,
    cores=1,
    ...
){
    
    #input and input checks
    counts <- getData(spCounts, "counts.cpm")
    
    #calculate fractions
    groupMeans <- getData(spUnsupervised, "groupMeans")
    fractions <- rep(1.0/(dim(groupMeans)[2]), (dim(groupMeans)[2]))
    
    #subset top genes for use with optimization
    selectInd <- getData(spUnsupervised, "selectInd")
    cellTypes <- groupMeans[selectInd, ]
    slice <- counts[selectInd, ]
    
    ##run pySwarm
    result <- .runPyRSwarm(
        cellTypes=cellTypes,
        slice=slice,
        fractions=fractions,
        distFun=distFun,
        maxiter=maxiter,
        swarmsize=swarmsize,
        cores=cores
    )
    
    #create object
    new("spSwarm",
        spSwarm=result,
        arguments = list(
            maxiter=maxiter,
            swarmsize=swarmsize
        )
    )
})


##run optimization
.runPyRSwarm <- function(
    cellTypes,
    slice,
    fractions,
    distFun,
    maxiter,
    swarmsize,
    cores
){
    
    control=list(maxit=maxiter, s=swarmsize)
    
    optim.fn <- function(i) {
        oneslice <- slice[,i]
        cat("running", i, "\n")
        psoptim(
            par=fractions,
            fn=distFun,
            cellTypes=cellTypes,
            slice=oneslice,
            lower=0,
            upper=1,
            control=control)$par
    }
    
    result <- as.data.frame(mclapply(1:(dim(slice)[2]), optim.fn, mc.cores=cores))
    output <- as.data.frame(t(result))
    colnames(output) <- colnames(cellTypes)
    rownames(output) <- colnames(slice)
    return(output)
}

.makeSyntheticSlice <- function(celltypes, fractions) {
    return(colSums(t(celltypes) * fractions))
}

# Various dist functions. Probably better to use match.arg and not export (so as to avoid cluttering the namespace), but leaving it like this for now.

#' @export
distToSlice <- function(fractions, cellTypes, slice) {
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum(abs(a - slice)) # Abs dist
}

#' @export
distToSliceEuclid <- function(fractions, cellTypes, slice) {
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum((a - slice)^2) # Euclidean dist
}

#' @export
distToSlicePearson <- function(fractions, cellTypes, slice) {
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum(1-(cor(a,slice))) # Pearson 'dist'
}

#' @export
exportByPoisson <- function(swarm, edge.cutoff, min.pval=0, min.num.edges=0) {
    swarmT <- as.data.frame(t(swarm@spRSwarm)) # FIXME: transpose because of cross-commit, should just fix indices below istead.
    sobj.bool <- swarmT > edge.cutoff
    nodes <- rownames(swarmT)
    nclusts <- length(nodes)
    edges <- data.frame(node1=rep(nodes, each=nclusts), node2=rep(nodes, nclusts), strength=rep(0, nclusts^2))
    for(i in 1:(dim(sobj.bool)[2])) {
        o <- which(sobj.bool[,i])
        if(length(o) > 1) {
            cat("o:", o, "\n")
            for(j in 1:(length(o)-1)) {
                for(k in (j+1):length(o)) {
                    cat(o[j], "\t", o[k], "\t")
                    ind1 <- (o[j]-1)*nclusts+o[k]
                    cat(ind1, "\t")
                    edges[ind1,3] <- edges[ind1,3]+1
                    ind2 <- (o[k]-1)*nclusts+o[j]
                    cat(ind2, "\n")
                    edges[ind2,3] <- edges[ind2,3]+1
                }
            }
        }
    }
    mean.edges <- mean(edges$strength)
    edges$pval <- ppois(edges$strength, mean.edges)
    edges <- edges[edges[[4]] >= min.pval & edges[[3]] >= min.num.edges,]
    edges
}
