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
#' @param spCounts an spCount object with multiplets.
#' @param spUnsupervised an spCount object.
#' @param distFun The distance function used to calculate the cost.
#' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of swarm particals.
#' @param cores The number of cores to be used while running spRSwarm.
#' @param spSwarm The spSwarm results.
#' @param costs The costs after optimization.
#' @param convergence The convergence output from psoptim. One value per multiplet.
#' @param arguments Arguments passed to the spRSwarm function.
#' @param object spRSwarm object.
#' @param n Data to extract from spRSwarm object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spSwarm output.
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
    seed=11,
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
    multiplets <- counts[selectInd, ]
    
    ##run pySwarm
    tmp <- .runPyRSwarm(
        cellTypes=cellTypes,
        multiplets=multiplets,
        fractions=fractions,
        distFun=distFun,
        maxiter=maxiter,
        swarmsize=swarmsize,
        cores=cores,
        seed=seed
    )
    result <- tmp[[1]]
    cost <- tmp[[2]]
    convergence <- tmp[[3]]
    
    #create object
    new("spSwarm",
        spSwarm=result,
        costs=cost,
        convergence=convergence,
        arguments = list(
            maxiter=maxiter,
            swarmsize=swarmsize
        )
    )
})


##run optimization
.runPyRSwarm <- function(
    cellTypes,
    multiplets,
    fractions,
    distFun,
    maxiter,
    swarmsize,
    cores,
    seed
){
    
    control=list(maxit=maxiter, s=swarmsize)
    
    set.seed(seed)
    tmp <- mclapply(
        1:(dim(multiplets)[2]),
        function(i)
            .optim.fn(
                i,
                fractions,
                distFun,
                cellTypes,
                control,
                multiplets
            ),
        mc.cores=cores
    )
    
    #compile results
    output <- data.frame(t(sapply(tmp, function(j) j[[1]])))
    cost <- sapply(tmp, function(j) j[[2]])
    counts <- t(sapply(tmp, function(j) j[[3]]))
    convergence <- sapply(tmp, function(j) j[[4]])
    convergenceKey <- c(
        "Maximal number of function evaluations reached." = 1,
        "Maximal number of iterations reached." = 2,
        "Maximal number of restarts reached." = 3,
        "Maximal number of iterations without improvement reached." = 4
    )
    convergence <- names(convergenceKey)[match(convergence, convergenceKey)]
    
    
    #normalize swarm output
    output <- output * 1/rowSums(output)

    colnames(output) <- colnames(cellTypes)
    rownames(output) <- colnames(multiplets)
    return(list(output, cost, convergence))
}

.optim.fn <- function(
    i,
    fractions,
    distFun,
    cellTypes,
    control,
    multiplets
){
    oneMultiplet <- multiplets[,i]
    psoptim(
        par=fractions,
        fn=distFun,
        cellTypes=cellTypes,
        oneMultiplet=oneMultiplet,
        lower=0,
        upper=1,
        control=control
    )
}

#calculates the sum for each gene after each cell type has been multiplied by
#fractions
.makeSyntheticSlice <- function(
    cellTypes,
    fractions
){
    return(colSums(t(cellTypes) * fractions))
}

# Various dist functions. Probably better to use match.arg and not export
#(so as to avoid cluttering the namespace), but leaving it like this for now.

#' @export
distToSlice <- function(
    fractions,
    cellTypes,
    oneMultiplet
){
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum(abs(a - oneMultiplet))
}

#' @export
distToSliceTopFour <- function(
    fractions,
    cellTypes,
    oneMultiplet
){
    if(sum(fractions) == 0) {
        return(999999999)
    }
    fraction[fraction < sort(fraction, decresing=T)[4]] <- 0 # 4 is an arbitrary number...
    cat(fraction)
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum(abs(a - oneMultiplet))
}

#' @export
distToSliceEuclid <- function(
    fractions,
    cellTypes,
    oneMultiplet
){
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum((a - oneMultiplet)^2)
}

#' @export
distToSlicePearson <- function(
    fractions,
    cellTypes,
    oneMultiplet
){
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    1-(cor(a, oneMultiplet))
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

.funx <- function(
    row,
    totcomb
){
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


#' calcResiduals
#'
#' Subtitle
#'
#' Calculates the residuals for each gene and multiplet after deconvolution
#' based on the spSwarm results.
#'
#' @name calcResiduals
#' @rdname calcResiduals
#' @aliases calcResiduals
#' @param spCounts An spCounts object with multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @param spSwarm An spSwarm object.
#' @param clusters A character vector of length 2 indicating 2 classes to
#'    specifically extract residuals from.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param object spRSwarm object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spSwarm connection strengths and p-values.
#' @author Jason T. Serviss
#' @keywords calcResiduals
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname calcResiduals
#' @export

calcResiduals <- function(
    spCounts,
    spUnsupervised,
    spSwarm,
    clusters=NULL,
    edge.cutoff=NULL,
    distFun = function(frac, multiplets) {sum(abs(frac - multiplets))},
    ...
){
    #spCounts should only include multiplets
    
    frac <- getData(spSwarm, "spSwarm")
    groupMeans <- getData(spUnsupervised, "groupMeans")
    selectInd <- getData(spUnsupervised, "selectInd")
    counts <- getData(spCounts, "counts.cpm")
    
    cellTypes <- groupMeans[selectInd, ]
    multiplets <- counts[selectInd, ]
    
    a <- sapply(1:nrow(frac), function(j) .makeSyntheticSlice(cellTypes, j))
    diff <- distFun(a, multiplets)
    #diff <- diff * 1/rowSums(diff)

    colnames(diff) <- colnames(counts)
    
    if(!is.null(clusters) & !is.null(edge.cutoff)) {
        diff <- diff[ , selectClustersOnEdge(
            spSwarm,
            edge.cutoff,
            clusters[1],
            clusters[2]
        )]
    }
    return(diff)
}

#.makeSyntheticSlice <- function(
#    celltypes,
#    fractions
#){
#    return(colSums(t(cellTypes) * fractions))
#}


#' selectClustersOnEdge
#'
#' Subtitle
#'
#' Description
#'
#' @name selectClustersOnEdge
#' @rdname selectClustersOnEdge
#' @aliases selectClustersOnEdge
#' @param spSwarm An spSwarm object.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param clust1 A character vector of length 1 indicating one node which forms
#'    an edge with clust2 for which multiplets contributing to that edge will be
#'    reported.
#' @param clust2 A character vector of length 1 indicating one node which forms
#'    an edge with clust1 for which multiplets contributing to that edge will be
#'    reported.
#' @param object spRSwarm object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spSwarm connection strengths and p-values.
#' @author Jason T. Serviss
#' @keywords selectClustersOnEdge
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname selectClustersOnEdge
#' @export

selectClustersOnEdge <- function(
    spSwarm,
    edge.cutoff,
    clust1,
    clust2
){
    frac <- getData(spSwarm, "spSwarm")[,c(clust1, clust2)]
    o <- apply(frac, 1, function(x) {all(x > edge.cutoff)})
    return(rownames(frac)[o])
}