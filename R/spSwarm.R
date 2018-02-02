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
#' @param cellNumbers Tibble; Output from \code{estimateCells} function.
#' @param e Numeric; The epsilon value for the .complexityPenilty unit.
#' @param spSwarm The spSwarm results.
#' @param costs The costs after optimization.
#' @param convergence The convergence output from psoptim. One value per
#'    multiplet.
#' @param stats The stats output from psoptim.
#' @param arguments Arguments passed to the spSwarm function.
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
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#' cn <- estimateCells(cObjSng, cObjMul)
#' sObj <- spSwarm(cObjMul, testUns, distFun = "dtsnCellNum", cellNumbers = cn, e = 0.0025)
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
  distFun = "distToSlice",
  maxiter = 10,
  swarmsize = 150,
  cores = 1,
  seed = 11,
  norm = TRUE,
  report = FALSE,
  reportRate = NULL,
  cellNumbers = NULL,
  e = NULL,
  ...
){
    
    #put a check here to make sure all slots in the spUnsupervised object are filled.
    
    if(length(distFun) != 1) {
        stop("Please provide a valid distFun argument.")
    }
    if(distFun == "dtsnCellNum" & (is.null(cellNumbers) | is.null(e))) {
      stop("cellNumbers and e must be provided with dtsnCellNum distFun.")
    }
    
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
        cellNumbers = cellNumbers,
        e = e,
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
    cellNumbers,
    e,
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
    
    if(!is.null(cellNumbers)) {
      matchIdx <- match(colnames(multiplets), cellNumbers$sampleName)
      cellNumbers <- cellNumbers$cellNumberMedian[matchIdx]
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
                cellNumbers,
                e,
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
    cellNumbers,
    e,
    ...
){
    oneMultiplet <- multiplets[, i]
    cellNumber <- cellNumbers[i]
    
    psoptim(
        par = fractions,
        fn = distFun,
        cellTypes = cellTypes,
        oneMultiplet = oneMultiplet,
        lower = 0,
        upper = 1,
        control = control,
        i = i,
        cellNumber = cellNumber,
        e = e,
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

#function which calculates the complexity penalty
.complexityPenilty <- function(k, e, cellNumber) {
  n <- k / log(cellNumber)
  u <- n * e
  1 + u
}

# Various dist functions. Probably better to use match.arg and not export
#(so as to avoid cluttering the namespace), but leaving it like this for now.

dtsnCellNum <- function(
    fractions,
    cellTypes,
    oneMultiplet,
    e,
    cellNumber,
    ...
){
  if(sum(fractions) == 0) {
      return(999999999)
  }
  normFractions <- fractions / sum(fractions)
  cellTypes <- cellTypes/mean(cellTypes)
  a <- .makeSyntheticSlice(cellTypes, normFractions)
  a <- a/mean(a)
  k <- length(which(normFractions > 0))
  penalty <- .complexityPenilty(k, e, cellNumber)
  sum(abs((oneMultiplet - a) / (a+1))) * penalty
}

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
    n <- length(fractions)
    e <- sum(abs(a - oneMultiplet)^2) * 1/n
    k <- length(which(fractions > 0))
    
    (n * log(e)) + (k * log(n))
}

bicLinear <- function(
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
    e <- sum(abs(a - oneMultiplet)) * 1/length(fractions)
    n <- length(fractions)
    k <- length(which(fractions > 0))
    
    (n * log(e)) + (k * log(n))
}

bic10 <- function(
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
    n <- length(fractions)
    e <- sum(abs(a - oneMultiplet)^2) * 1/n
    k <- length(which(fractions > 0))
    
    (n * log(e)) + (k * log10(n))
}

bic10.2 <- function(
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
    n <- length(fractions)
    e <- sum(abs(a - oneMultiplet)^2) * 1/n
    k <- length(which(fractions > 0))
    
    (n * log(e)) + ((k / n) * log10(n))
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
#' output <- spSwarmPoisson(testSwa, 1/10.5)
#'
#'
NULL

#' @rdname spSwarmPoisson
#' @importFrom stats ppois
#' @importFrom dplyr filter pull rowwise do ungroup
#' @importFrom rlang .data
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
        mean <- edges %>%
            filter(.data$from %in% c(f, t) | .data$to %in% c(f, t)) %>%
            pull(weight) %>%
            mean()
        
        ppois(weight, mean, lower.tail = FALSE) %>%
            as_tibble() %>%
            rename(pval = .data$value)
    }
    
    edges %>%
        rowwise() %>%
        do(
            bind_cols(
                .data,
                ps(edges, .data$from, .data$to, .data$weight)
            )
        ) %>%
        ungroup %>%
        filter(.data$pval <= min.pval & .data$weight >= min.num.edges)
}

.calculateWeight <- function(
    mat,
    logic,
    ...
){
  homo <- paste(sort(colnames(mat)), sort(colnames(mat)), sep = "-")
  hetero <- apply(
    apply(combn(colnames(mat), 2), 2, sort),
    2,
    paste,
    collapse = "-"
  )
  totcomb <- c(homo, hetero)
    
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
          apply(combn(pick, 2), 2, sort),
          2,
          paste,
          collapse = "-"
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
#' s <- grepl("^s", colnames(testCounts))
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#' output <- calcResiduals(cObjMul, testUns, testSwa)
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
#' output <- getMultipletsForEdge(testSwa, 1/10.5, data.frame("A1", "B1"))
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
#' @importFrom rlang .data
#' @importFrom dplyr mutate select rename pull
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
        select(.data$variables, .data$from, .data$to, -.data$names) %>%
        rename(multiplet = .data$variables)

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
#' output <- getEdgesForMultiplet(testSwa, 1/10.5, "m.A1B1")
#'
NULL

#' @rdname getEdgesForMultiplet
#' @export
#' @importFrom rlang .data
#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select
#' @importFrom stats setNames

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
    frac <- getData(spSwarm, "spSwarm")[multiplet, ]
    map_dfr(multiplet, .edgeFunx, edge.cutoff, frac, s)
})

.edgeFunx <- function(x, edge.cutoff, frac, s) {
    keep <- as.logical(frac[x, ] > edge.cutoff)
    
    if(length(which(keep)) > 1) {
        
        combs <- combn(names(frac)[keep], 2)
        sapply(1:ncol(combs), function(y) {
            bool1 <- s$from == combs[1, y] & s$to == combs[2, y]
            bool2 <- s$from == combs[2, y] & s$to == combs[1, y]
            as.character(distinct(filter(s, bool1 | bool2)[, 1:2]))
        }) %>%
        t() %>%
        as_tibble() %>%
        setNames(c("from", "to")) %>%
        add_column("multiplet" = x) %>%
        select(.data$multiplet, .data$from, .data$to)
        
    } else if(length(which(keep)) == 1) {
        n <- names(frac)[keep]
        filter(s, .data$to == n & .data$from == n) %>%
        add_column("multiplet" = x) %>%
        select(-.data$weight, -.data$pval)
    }
}
