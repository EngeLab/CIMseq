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
#' @param seed The desired seed to set before running.
#' @param norm Logical indicating if the sum of fractions should equal 1.
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

setGeneric("spSwarm", function(
    spCounts,
    spUnsupervised,
    ...
){
    standardGeneric("spSwarm")
})


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
    norm=TRUE,
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
    multiplets <- matrix(
        counts[selectInd, ],
        ncol=ncol(counts),
        dimnames=list(1:length(selectInd), colnames(counts))
    )
    
    ##run pySwarm
    tmp <- .runPyRSwarm(
        cellTypes = cellTypes,
        multiplets = multiplets,
        fractions = fractions,
        distFun = distFun,
        maxiter = maxiter,
        swarmsize = swarmsize,
        cores = cores,
        seed = seed,
        norm = norm,
        ...
    )
    result <- tmp[[1]]
    cost <- tmp[[2]]
    convergence <- tmp[[3]]
    #stats <- tmp[[4]]
    
    #create object
    new("spSwarm",
        spSwarm=result,
        costs=cost,
        convergence = convergence,
        #stats = stats,
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
    seed,
    norm,
    ...
){
    
    control=list(
        maxit=maxiter,
        s=swarmsize
        #trace=1,
        #REPORT=10,
        #trace.stats=TRUE
    )
    
    set.seed(seed)
    to <- if(ncol(multiplets) == 1) {to <- 1} else {to <- dim(multiplets)[2]}
    
    tmp <- mclapply(
        1:to,
        function(i)
            .optim.fn(
                i,
                fractions,
                distFun,
                cellTypes,
                control,
                multiplets,
                ...
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
    #stats <- lapply(tmp, function(j) j[[6]])
    
    #normalize swarm output
    if(norm) {
        output <- output * 1/rowSums(output)
    }
    
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
    multiplets,
    ...
){
    oneMultiplet <- multiplets[,i]
    psoptim(
        par=fractions,
        fn=distFun,
        cellTypes=cellTypes,
        oneMultiplet=oneMultiplet,
        lower=0,
        upper=1,
        control=control,
        i=i,
        ...
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
    oneMultiplet,
    ...
){
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    sum(abs(a - oneMultiplet))
}

#' @export
distToSliceNorm <- function(
    fractions,
    cellTypes,
    oneMultiplet,
    i,
    ...
){
    if(sum(fractions) == 0) {
        return(999999999)
    }
    normFractions <- fractions / sum(fractions)
    cellTypes <- cellTypes/mean(cellTypes)
    a = .makeSyntheticSlice(cellTypes, normFractions)
    a <- a/mean(a)
    sum(abs((oneMultiplet - a) / (a+1)))
}

#' @export
distToSliceTop <- function(
    fractions,
    cellTypes,
    oneMultiplet,
    i,
    ...
){
    if("cells" %in% names(list(...))) {
        l <- list(...)
        cells <- ceiling(l[['cells']][i])
    }
    cat(fractions, "\n")
    if(sum(fractions) == 0) {
        return(999999999)
    }
    fractions[fractions < sort(fractions, decreasing=TRUE)[cells]] <- 0
    a = .makeSyntheticSlice(cellTypes, fractions)
    sum(abs(a - oneMultiplet))
}

#' @export
distToSliceEuclid <- function(
    fractions,
    cellTypes,
    oneMultiplet,
    i,
    ...
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
    oneMultiplet,
    i,
    ...
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
#' @importFrom stats ppois
#' @export

spSwarmPoisson <- function(
    spSwarm,
    edge.cutoff,
    min.pval=1,
    min.num.edges=0,
    ...
){
    mat <- getData(spSwarm, "spSwarm")
    logic <- .fractionCutoff(mat, edge.cutoff)
    
    #calcluate weight
    edges <- .calculateWeight(
        mat,
        logic
    )
    
    #calculate p-value
    out <- .calculateP(
        edges,
        min.pval,
        min.num.edges
    )
    
    return(out)
}

.fractionCutoff <- function(mat, cutoff) {
    if(length(cutoff) == 1) {
        return(mat > cutoff)
    } else {
        logic <- t(sapply(1:nrow(mat), function(j) {
            mat[j,] >= as.numeric(
                sort(mat[j,], decreasing=TRUE)[ceiling(cutoff[j])]
            )
        }))
        colnames(logic) <- colnames(mat)
        return(logic)
    }
}
.calculateP <- function(
    edges,
    min.pval,
    min.num.edges,
    ...
){
    means <- sapply(1:nrow(edges), function(o) {
        mean(subset(
            edges,
            from %in% c(edges[o, 1], edges[o, 2]) |
            to %in% c(edges[o, 1], edges[o, 2])
        )$weight)
    })
    
    edges$pval <- ppois(edges$weight, means, lower.tail=FALSE)
    out <- edges[edges$pval < min.pval & edges$weight >= min.num.edges, ]
    return(out)
}

.calculateWeight <- function(
    mat,
    logic,
    ...
){
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
    
    edges <- data.frame(
        from=xy[,1],
        to=xy[,2],
        weight=res,
        stringsAsFactors=FALSE,
        row.names=1:nrow(xy)
    )
    
    return(edges)
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
#' @param distFun The distance function to be used to calculate the residuals.
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
    distFun = function(frac, multiplets){(abs(multiplets - frac) / (frac + 1))},
    ...
){
    #spCounts should only include multiplets
    
    frac <- getData(spSwarm, "spSwarm")
    groupMeans <- getData(spUnsupervised, "groupMeans")
    selectInd <- getData(spUnsupervised, "selectInd")
    counts <- getData(spCounts, "counts.cpm")
    
    cellTypes <- groupMeans[selectInd, ]
    multiplets <- counts[selectInd, ]
    multiplets <- multiplets/mean(multiplets)
    
    a <- sapply(1:nrow(frac), function(j)
        .makeSyntheticSlice(cellTypes, as.numeric(frac[j,]))
    )
    colnames(a) <- rownames(frac)
    a <- a/mean(a)
    
    diff <- sapply(1:ncol(a), function(x)
        distFun(a[,x],multiplets[,x])
    )

    colnames(diff) <- rownames(frac)
    
    if(!is.null(clusters) & !is.null(edge.cutoff)) {
        diff <- diff[ , getMultipletsForEdge(
            spSwarm,
            edge.cutoff,
            clusters[1],
            clusters[2]
        )]
    }
    return(diff)
}

#' getMultipletsForEdge
#'
#' Returns the names of the multiplets that are associated with an edge.
#'
#' Description
#'
#' @name getMultipletsForEdge
#' @rdname getMultipletsForEdge
#' @aliases getMultipletsForEdge
#' @param spSwarm An spSwarm object.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param edges A data frame indicating the edges of interest. Edges are
#'    indicated by the nodes they connect with one node in column one and the
#'    other node in column 2.
#' @param ... additional arguments to pass on
#' @return If the edges argument only includes on row, a vector of multiplet
#'    names is returned. If several edges are interogated a list is returned
#'    with one element per edge containing the names of the multiplets.
#' @author Jason T. Serviss
#' @keywords getMultipletsForEdge
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname getMultipletsForEdge
#' @export

setGeneric("getMultipletsForEdge", function(
    spSwarm,
    ...
){
    standardGeneric("getMultipletsForEdge")
})

#' @rdname getMultipletsForEdge
#' @export

setMethod("getMultipletsForEdge", "spSwarm", function(
    spSwarm,
    edge.cutoff,
    edges,
    ...
){
    edges[,1] <- as.character(edges[,1])
    edges[,2] <- as.character(edges[,2])
    
    mulForEdges <- sapply(1:nrow(edges), function(j) {
        frac <- getData(spSwarm, "spSwarm")[,c(edges[j,1], edges[j,2])]
        o <- apply(frac, 1, function(x) {all(x > edge.cutoff)})
        rownames(frac)[o]
    })
    
    if(class(mulForEdges) == "matrix") {
        return(unlist(mulForEdges))
    } else {
        names(mulForEdges) <- paste(edges[,1], edges[,2], sep="-")
        return(mulForEdges)
    }
})

#' getEdgesForMultiplet
#'
#' Returns the names of the edges are associated with a multiplet.
#'
#' Description
#'
#' @name getEdgesForMultiplet
#' @rdname getEdgesForMultiplet
#' @aliases getEdgesForMultiplet
#' @param spSwarm An spSwarm object.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param multiplet The name of the multiplet of interest.
#' @param ... additional arguments to pass on
#' @return Edge names.
#' @author Jason T. Serviss
#' @keywords getEdgesForMultiplet
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname getEdgesForMultiplet
#' @export

setGeneric("getEdgesForMultiplet", function(
    spSwarm,
    ...
){
    standardGeneric("getEdgesForMultiplet")
})

#' @rdname getEdgesForMultiplet
#' @export

setMethod("getEdgesForMultiplet", "spSwarm", function(
    spSwarm,
    edge.cutoff,
    multiplet,
    ...
){
    frac <- getData(spSwarm, "spSwarm")[multiplet,]
    combs <- combn(names(frac)[frac > edge.cutoff], 2)
    s <- spSwarmPoisson(spSwarm, edge.cutoff=edge.cutoff)
    
    #out <- as.data.frame(t(sapply(
    #    1:ncol(combs),
    #    function(j)
    #        s[s$from == combs[1,j] & s$to == combs[2,j], ]
    #)))
    out <- s[s[,1] %in% t(combs)[,1] & s[,2] %in% t(combs)[,2], ]
    
    if(nrow(out) == ncol(combs)) {
        return(out)
    } else {
        stop("somethings went wrong, check the code")
    }
})

#' permuteSwarm
#'
#' Write something
#'
#' Description
#'
#' @name permuteSwarm
#' @rdname permuteSwarm
#' @aliases permuteSwarm
#' @param spSwarm An spSwarm object.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param multiplet The name of the multiplet of interest.
#' @param ... additional arguments to pass on
#' @return Edge names.
#' @author Jason T. Serviss
#' @keywords permuteSwarm
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname permuteSwarm
#' @export

setGeneric("permuteSwarm", function(
    spCounts,
    ...
){
    standardGeneric("permuteSwarm")
})

#' @rdname permuteSwarm
#' @export

setMethod("permuteSwarm", "spCounts", function(
    spCounts,
    spCountsMul,
    spUnsupervised,
    spSwarm,
    distFun = distToSlice,
    maxiter = 10,
    swarmsize = 150,
    cores=1,
    seed=11,
    norm=TRUE,
    iter,
    ...
){
    classes <- getData(spUnsupervised, "classification")
    
    permMatrix <- .makePermutations(classes, iter)
    
    rows <- ncol(getData(spCountsMul, "counts"))
    cols <- length(unique(classes))
    dims <- iter
    permData <- array(
        NA,
        dim=c(rows, cols, dims),
        dimnames=list(
            colnames(getData(spCountsMul, "counts")),
            sort(unique(classes)),
            1:iter
        )
    )
    
    for(i in 1:iter) {
        classification(spUnsupervised) <- permMatrix[,i]
        groupMeans(spUnsupervised) <- averageGroupExpression(
            spCounts,
            permMatrix[,i],
            weighted=FALSE
        )
        sObj <- spSwarm(
            spCountsMul,
            spUnsupervised,
            distFun = distToSlice,
            maxiter = maxiter,
            swarmsize = swarmsize,
            cores=cores,
            seed=seed,
            norm=norm
        )
        mat <- getData(sObj, "spSwarm")
        order <- as.matrix(mat[ , order(colnames(mat))])
        permData[,,i] <- order
    }
    
    return(.calculatePermP(permData, iter))
})

.makePermutations <- function(classes, iter){
    sapply(1:iter, function(j) {
        classes[sample(1:length(classes), size=length(classes), replace=FALSE)]
    })
}

.calculatePermP <- function(permData, spSwarm) {
    swarm <- getData(spSwarm, "spSwarm")
    logic <- .fractionCutoff(swarm, 0)
    edges <- .calculateWeight(swarm, logic, 0, 1, 0)[,1:3]
    
    #remove self connections for now
    edges <- edges[edges$from != edges$to, ]
    
    p <- c()
    #get null for each edge
    for(i in 1:nrow(edges)) {
        node1 <- edges[i, 1]
        node2 <- edges[i, 2]
        idx1 <- which(colnames(swarm) == node1)
        idx2 <- which(colnames(swarm) == node2)

        swarm1 <- swarm[,idx1]
        swarm2 <- swarm[,idx2]
        
        null1 <- permData[,idx1,]
        null2 <- permData[,idx2,]
        
        l1 <- swarm1 > null1
        l2 <- swarm2 > null2
        l <- l1+l2
        l[l < 2] <- 0
        l[l == 2] <- 1
        
        p <- c(p, sum(l)/(iter*nrow(swarm)))

    }
    
    edges$p <- p
    return(edges)
}

