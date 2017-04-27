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
#' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of swarm particals.
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
    tmp <- .runPyRSwarm(
        cellTypes=cellTypes,
        slice=slice,
        fractions=fractions,
        distFun=distFun,
        maxiter=maxiter,
        swarmsize=swarmsize,
        cores=cores
    )
    result <- tmp[[1]]
    cost <- tmp[[2]]
    #create object
    new("spSwarm",
        spSwarm=result,
        costSum=sum(cost),
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
    
    optim.fn <- function(i, fractions, distFun, cellTypes, control) {
        oneslice <- slice[,i]
        psoptim(
            par=fractions,
            fn=distFun,
            cellTypes=cellTypes,
            slice=oneslice,
            lower=0,
            upper=1,
            control=control)
    }
    
    
    tmp <- mclapply(
        1:(dim(slice)[2]),
        function(i)
            optim.fn(
                i,
                fractions,
                distFun,
                cellTypes,
                control
            ),
        mc.cores=cores
    )
    
    #compile results
    output <- data.frame()
    cost <- c()
    counts <- tmp[[1]][[3]]
    convergence <-ifelse(
        tmp[[1]][[4]] == 1,
        "Maximal number of function evaluations reached.",
        ifelse(
            tmp[[1]][[4]] == 2,
            "Maximal number of iterations reached.",
            ifelse(
                tmp[[1]][[4]] == 3,
                "Maximal number of restarts reached.",
                "Maximal number of iterations without improvement reached."
            )
        )
    )
    
    for(i in 1:length(tmp)) {
        curr <- tmp[[i]]
        output <- rbind(output, as.data.frame(t(curr[[1]])))
        cost <- c(cost, curr[[2]])
    }
    
    #normalize swarm output
    rs <- rowSums(output)
    output <-  t(apply(output, 1, function(x) x/rs[x]))

    colnames(output) <- colnames(cellTypes)
    rownames(output) <- colnames(slice)
    return(list(output, cost))
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
distToSliceTopFour <- function(fractions, cellTypes, slice) {
    if(sum(fractions) == 0) {
        return(999999999)
    }
    fraction[fraction < sort(fraction, decresing=T)[4]] <- 0 # 4 is an arbitrary number...
    cat(fraction)
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

#' spSwarmPoisson
#'
#' Subtitle
#'
#' Description
#'
#' @name spSwarmPoisson
#' @rdname spSwarmPoisson
#' @aliases spSwarmPoisson
#' @param spSwarm An spSwarm object.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param min.pval Minimum p-value to report.
#' @param min.num.edges Minimum number of observed edges to report a connection.
#' @param object spRSwarm object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spSwarm connection strengths and p-values.
#' @author Jason T. Serviss
#' @keywords spSwarmPoisson
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname spSwarmPoisson
#' @export

spSwarmPoisson <- function(
    spSwarm,
    edge.cutoff,
    min.pval=1,
    min.num.edges=0,
    ...
){
    mat <- getData(spSwarm, "spSwarm")
    logic <- mat > edge.cutoff
    
    totcomb <- c(
        paste(colnames(mat), colnames(mat), sep=""),
        apply(
            t(combn(colnames(mat),2)),
            1,
            paste,
            collapse=""
        )
    )
    
    com <- apply(logic, 1, .funx, totcomb)
    rownames(com) <- totcomb
    res <- apply(com, 1, sum)
    xy <- rbind(
        matrix(c(colnames(mat), colnames(mat)), ncol=2),
        t(combn(colnames(mat), 2))
    )
    edges <- data.frame(from=xy[,1], to=xy[,2], weight=res)
    
    mean.edges <- mean(edges$weight)
    edges$pval <- ppois(edges$weight, mean.edges, lower.tail=FALSE)
    out <- edges[edges$pval < min.pval & edges$weight >= min.num.edges, ]
    return(out)
}

.funx <- function(row, totcomb){
    pick <- names(row)[which(row)]
    if(length(pick) == 1) {
        out <- totcomb %in% paste(pick, pick, sep="")
    } else {
        out <- totcomb %in% apply(
            t(combn(pick, 2)),
            1,
            paste,collapse=""
        )
    }
    return(out)
}

#spSwarmPoisson <- function(
#    spSwarm,
#    edge.cutoff,
#    min.pval=1,
#    min.num.edges=0,
#    ...
#){
# FIXME: transpose because of cross-commit, should just fix indices below istead
#    swarmT <- as.data.frame(t(getData(spSwarm, "spSwarm")))
#    sobj.bool <- swarmT > edge.cutoff
#    nodes <- rownames(swarmT)
#    nclusts <- length(nodes)
#    edges <- data.frame(
#        from=rep(nodes, each=nclusts),
#        to=rep(nodes, nclusts),
#        weight=rep(0, nclusts^2)
#    )
    
    #calculate edge weight, i.e. number of observed edges
#    for(i in 1:(dim(sobj.bool)[2])) {
#        o <- which(sobj.bool[,i])
#        if(length(o) > 1) {
#            for(j in 1:(length(o)-1)) {
#                for(k in (j+1):length(o)) {
#                    ind1 <- (o[j]-1)*nclusts+o[k]
#                    edges[ind1,3] <- edges[ind1,3]+1
#                    ind2 <- (o[k]-1)*nclusts+o[j]
#                    edges[ind2,3] <- edges[ind2,3]+1
#                }
#            }
#        } else { #add self-connections
#            edges[o,3] <- edges[ind1,3]
#        }
#    }
    
    #calculate p-value
#    mean.edges <- mean(edges$weight)
#    edges$pval <- ppois(edges$weight, mean.edges, lower.tail=FALSE)
#    out <- edges[edges$pval < min.pval & edges$weight >= min.num.edges, ]
#    return(out)
#}
