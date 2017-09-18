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
#' @param distFun The distance function used to calculate the cost. Either the
#'    name of a custom function in the local environment or one of the included
#'    functions, i.e. \code{distToSlice, distToSliceNorm, distToSliceTop,
#'    distToSliceEuclid, distToSlicePearson, bic}.
#' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of swarm particals.
#' @param cores The number of cores to be used while running spRSwarm.
#' @param seed The desired seed to set before running.
#' @param norm Logical indicating if the sum of fractions should equal 1.
#' @param report Logical indicating if additional reporting from the
#'   optimization should be included.
#' @param reportRate If report is TRUE, the iteration interval that a report
#'    should be generated.
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
    distFun = c(
        "distToSlice",
        "distToSliceNorm",
        "distToSliceTop",
        "distToSliceEuclid",
        "distToSlicePearson",
        "bic"
    ),
    maxiter = 10,
    swarmsize = 150,
    cores = 1,
    seed = 11,
    norm = TRUE,
    report = FALSE,
    reportRate = NULL,
    ...
){
    
    distFun <- match.fun(distFun)
    
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
        ncol = ncol(counts),
        dimnames = list(1:length(selectInd), colnames(counts))
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
        report = report,
        reportRate = reportRate,
        ...
    )
    result <- tmp[[1]]
    cost <- tmp[[2]]
    convergence <- tmp[[3]]
    stats <- tmp[[4]]
    
    #create object
    new("spSwarm",
        spSwarm = result,
        costs = cost,
        convergence = convergence,
        stats = stats,
        arguments = list(
            maxiter = maxiter,
            swarmsize = swarmsize
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
    report,
    reportRate,
    ...
){
    if(report) {
        control <- list(
            maxit = maxiter,
            s = swarmsize,
            trace = 1,
            REPORT = reportRate,
            trace.stats = TRUE
        )
    } else {
        control <- list(
            maxit = maxiter,
            s = swarmsize
        )
        stats <- list()
    }
    
    
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
        mc.cores = cores
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
    
    if(report) stats <- lapply(tmp, function(x) x[[6]])
    
    
    #normalize swarm output
    if(norm) {
        output <- output * 1/rowSums(output)
    }
    
    colnames(output) <- colnames(cellTypes)
    rownames(output) <- colnames(multiplets)
    return(list(output, cost, convergence, stats))
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
        par = fractions,
        fn = distFun,
        cellTypes = cellTypes,
        oneMultiplet = oneMultiplet,
        lower = 0,
        upper = 1,
        control = control,
        i = i,
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

bic <- function(
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
    a <- .makeSyntheticSlice(cellTypes, normFractions)
    e <- sum(abs(a - oneMultiplet)^2) * 1/length(fractions)
    n <- length(fractions)
    k <- length(which(fractions > 0))
    
    (n * log(e)) + (k * log(n))
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
#' @importFrom dplyr filter pull rowwise do ungroup
#' @export

spSwarmPoisson <- function(
    spSwarm,
    edge.cutoff,
    min.pval = 1,
    min.num.edges = 0,
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
    ps <- function(edges, f, t, weight) {
        edges %>%
            filter(from %in% c(f, t) | to %in% c(f, t)) %>%
            pull(weight) %>%
            mean() %>%
            ppois(weight, ., lower.tail = FALSE)
    }
    
    out <- edges %>%
        rowwise() %>%
        do(bind_cols(., tibble(pval = ps(edges, .$from, .$to, .$weight)))) %>%
        ungroup %>%
        filter(pval < min.pval & weight >= min.num.edges)
        
    #means <- sapply(1:nrow(edges), function(o) {
    #    mean(subset(
    #        edges,
    #        from %in% c(edges[o, 1], edges[o, 2]) |
    #        to   %in% c(edges[o, 1], edges[o, 2])
    #    )$weight)
    #})
    
    #edges$pval <- ppois(edges$weight, means, lower.tail=FALSE)
    #out <- edges[edges$pval < min.pval & edges$weight >= min.num.edges, ]
    return(out)
}

.calculateWeight <- function(
    mat,
    logic,
    ...
){
    totcomb <- c(
        paste(sort(colnames(mat)), sort(colnames(mat)), sep = ""),
        apply(
            t(combn(sort(colnames(mat)), 2)),
            1,
            paste,
            collapse = ""
        )
    )
    
    com <- apply(logic, 1, .funx, totcomb)
    rownames(com) <- totcomb
    res <- apply(com, 1, sum)
    xy <- rbind(
        matrix(c(sort(colnames(mat)), sort(colnames(mat))), ncol = 2),
        t(combn(sort(colnames(mat)), 2))
    )
    
    edges <- data.frame(
        from = xy[,1],
        to = xy[,2],
        weight = res,
        stringsAsFactors = FALSE,
        row.names = 1:nrow(xy)
    )
    
    return(edges)
}

.funx <- function(
    row,
    totcomb
){
    pick <- names(row)[which(row)]
    if(length(pick) == 1) {
        out <- totcomb %in% paste(pick, pick, sep = "")
    } else {
        out <- totcomb %in% apply(
            t(combn(pick, 2)),
            1,
            paste,collapse = ""
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
    clusters = NULL,
    edge.cutoff = NULL,
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
    edges[,1] <- as.character(pull(edges, 1))
    edges[,2] <- as.character(pull(edges, 2))
    
    mulForEdges <- lapply(1:nrow(edges), function(j) {
        cols <- c(pull(edges, 1)[j], pull(edges, 2)[j])
        
        if(identical(cols[1], cols[2])) {
            .self(j, cols, edges, spSwarm, edge.cutoff)
        } else {
            .nonSelf(j, cols, edges, spSwarm, edge.cutoff)
        }
        
    })
    
    names(mulForEdges) <- paste(pull(edges, 1), pull(edges, 2), sep = "-")
    
    namedListToTibble(mulForEdges) %>%
        mutate(from = gsub("(.*)-.*", "\\1", names)) %>%
        mutate(to = gsub(".*-(.*)", "\\1", names)) %>%
        select(variables, from, to, -names) %>%
        rename(multiplet = variables)

})

#note to self: when considering self-connections in multiplets it is important
#to remember that, at least now, if another edge is detected in the multiplet,
#the self-coonection will not be detected by the swarm optimization. For
#example, if the multiplet contains cell types A1, A1 and B1, the A1 self-
#connection will not be detected. On the other hand if the multiplet contains
#A1, A1, A1, then a self-connection will be reported.

.self <- function(j, cols, edges, spSwarm, edge.cutoff) {
    mat <- getData(spSwarm, "spSwarm")
    logic <- mat > edge.cutoff
    rs <- rowSums(logic)
    
    frac <- getData(spSwarm, "spSwarm")[, unique(cols)]
    o <- which(frac > edge.cutoff & rs == 1)
    rownames(getData(spSwarm, "spSwarm"))[o]
}

.nonSelf <- function(j, cols, edges, spSwarm, edge.cutoff) {
    frac <- getData(spSwarm, "spSwarm")[, cols]
    o <- apply(frac, 1, function(x) {all(x > edge.cutoff)})
    rownames(frac)[o]
}

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
#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select

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
    s <- spSwarmPoisson(spSwarm, edge.cutoff = edge.cutoff)
    frac <- getData(spSwarm, "spSwarm")[multiplet,]
    map_dfr(multiplet, .edgeFunx, edge.cutoff, frac, s)
})

.edgeFunx <- function(x, edge.cutoff, frac, s) {
    keep <- as.logical(frac[x, ] > edge.cutoff)
    
    if(length(which(keep)) > 1) {
        
        combs <- combn(names(frac)[keep], 2)
        sapply(1:ncol(combs), function(y) {
            bool1 <- s$from == combs[1, y] & s$to == combs[2, y]
            bool2 <- s$from == combs[2, y] & s$to == combs[1, y]
            as.character(filter(s, bool1 | bool2)[, 1:2])
        }) %>%
        t() %>%
        as_tibble() %>%
        setNames(c("from", "to")) %>%
        add_column(multiplet = x) %>%
        select(multiplet, from, to)
        
    } else if(length(which(keep)) == 1) {
        names(frac)[keep]
        filter(s, to == names(frac)[keep] & from == names(frac)[keep]) %>%
        add_column(multiplet = x) %>%
        select(-weight, -pval)
    }
}

#' permuteSwarm
#'
#' Write something
#'
#' Description
#'
#' @name permuteSwarm
#' @rdname permuteSwarm
#' @aliases permuteSwarm
#' @param spCountsSng An spCounts object containing singlets.
#' @param spCountsMul An spCounts object containing multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @param spSwarm An spSwarm object.
#' @param distFun The distance function used to calculate the cost. Either the
#'    name of a custom function in the local environment or one of the included
#'    functions, i.e. \code{distToSlice, distToSliceNorm, distToSliceTop,
#'    distToSliceEuclid, distToSlicePearson, bic}.
#' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of swarm particals.
#' @param cores The number of cores to be used while running spRSwarm.
#' @param seed The desired seed to set before running.
#' @param norm Logical indicating if the sum of fractions should equal 1.
#' @param iter The number of permutations to perform.
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
    spCountsSng,
    ...
){
    standardGeneric("permuteSwarm")
})

#' @rdname permuteSwarm
#' @export

setMethod("permuteSwarm", "spCounts", function(
    spCountsSng,
    spCountsMul,
    spUnsupervised,
    spSwarm,
    distFun = c(
        "distToSlice",
        "distToSliceNorm",
        "distToSliceTop",
        "distToSliceEuclid",
        "distToSlicePearson",
        "bic"
    ),
    maxiter = 10,
    swarmsize = 150,
    cores = 1,
    seed = 11,
    norm = TRUE,
    iter,
    ...
){
    distFun <- match.fun(distFun)

    classes <- getData(spUnsupervised, "classification")
    
    permMatrix <- .makePermutations(classes, iter, seed)
    permData <- .runPermutations(
        permMatrix,
        iter,
        spCounts,
        spCountsMul,
        spUnsupervised,
        classes,
        distFun = distToSlice,
        maxiter = maxiter,
        swarmsize = swarmsize,
        cores = cores,
        seed = seed,
        norm = norm
    )
    
    return(list(.calculatePermP(permData, spSwarm, iter), permData))
})

.makePermutations <- function(classes, iter, seed){
    set.seed(seed)
    sapply(1:iter, function(j) {
        classes[sample(1:length(classes), size=length(classes), replace=FALSE)]
    })
}

.runPermutations <- function(
    permMatrix,
    iter,
    spCountsSng,
    spCountsMul,
    spUnsupervised,
    classes,
    distFun = distToSlice,
    maxiter = 10,
    swarmsize = 150,
    cores = 1,
    seed = 11,
    norm = TRUE
){
    
    permData <- list()
    
    for(i in 1:iter) {
        classification(spUnsupervised) <- permMatrix[,i]
        groupMeans(spUnsupervised) <- averageGroupExpression(
            spCountsSng,
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
        permData[[i]] <- sObj
    }
    
    pEdges <- lapply(1:length(permData), function(j) {
        mat <- getData(permData[[j]], "spSwarm")
        logic <- .fractionCutoff(mat, 0)
        .calculateWeight(mat, logic)$weight
    })
    
    pEdges <- t(do.call(rbind,pEdges))
    return(pEdges)
}

.calculatePermP <- function(pEdges, spSwarm, iter) {
    rEdges <- spSwarmPoisson(spSwarm, edge.cutoff=0)
    rEdges$pval <- NULL
    
    rownames(pEdges) <- paste(rEdges$from, rEdges$to, sep="-")
    
    pValues <- sapply(1:nrow(rEdges), function(g) {
        if(rEdges[g, "weight"] == 0) {
            NA
        } else if(sum(pEdges[g,] >= rEdges[g, "weight"]) == 0) {
            10^-log10(iter)
        } else {
            sum(pEdges[g,] >= rEdges[g, "weight"]) / iter
        }
    })
    
    rEdges$pValue <- pValues
    return(rEdges)
}

#rows <- ncol(getData(spCountsMul, "counts"))
#cols <- length(unique(classes))
#dims <- iter
#permData <- array(
#    NA,
#    dim=c(rows, cols, dims),
#    dimnames=list(
#        colnames(getData(spCountsMul, "counts")),
#        sort(unique(classes)),
#        1:iter
#    )
#)

#mat <- getData(sObj, "spSwarm")
#order <- as.matrix(mat[ , order(colnames(mat))])
#permData[,,i] <- order



#.calculatePermP <- function(permData, spSwarm, edgeAlpha) {
#    swarm <- getData(spSwarm, "spSwarm")
#    logic <- .fractionCutoff(swarm, 0)
#    edges <- .calculateWeight(swarm, logic, 0, 1, 0)[,1:2]
#
#    #remove self connections for now
#    edges <- edges[edges$from != edges$to, ]
#
#    p <- c()
#    weight <- c()
#
#    #get null for each edge
#    for(i in 1:nrow(edges)) {
#        node1 <- edges[i, 1]
#        node2 <- edges[i, 2]
#        idx1 <- which(colnames(swarm) == node1)
#        idx2 <- which(colnames(swarm) == node2)
#
#        swarm1 <- swarm[,idx1]
#        swarm2 <- swarm[,idx2]
#
#        null1 <- permData[,idx1,]
#        null2 <- permData[,idx2,]
#
#        l1 <- swarm1 > null1
#        l2 <- swarm2 > null2
#        l <- l1+l2
#        l[l < 2] <- 0
#        l[l == 2] <- 1
#
#        p <- c(p, sum(l)/(iter*nrow(swarm)))
#
#        #calculate edge weight
#        l1 <- swarm1 < null1
#        l2 <- swarm2 < null2
#
#        n1 <- rowSums(l1)
#        n2 <- rowSums(l2)
#        n1P <- ifelse(n1 == 0, 10^-log10(iter), n1/iter)
#        n2P <- ifelse(n2 == 0, 10^-log10(iter), n2/iter)
#        weight <- c(weight, length(
#        intersect(
#        which(n1P < edgeAlpha),
#        which(n2P < edgeAlpha)
#        )
#        ))
#    }
#
#    edges$p <- p
#    edges$weight <- weight
#    return(edges)
#}




