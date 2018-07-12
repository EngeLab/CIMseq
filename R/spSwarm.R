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
#' @param spCountsSng an spCount object with singlets.
#' @param spCountsMul an spCount object with multiplets.
#' @param spUnsupervised an spCount object.
#' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of swarm particals.
#' @param nSyntheticMultiplets Numeric value indicating the number of synthetic
#'  multiplets to generate during deconvolution.
#' @param seed The desired seed to set before running.
#' @param norm Logical indicating if the sum of fractions should equal 1.
#' @param report Logical indicating if additional reporting from the
#'   optimization should be included.
#' @param reportRate If report is TRUE, the iteration interval that a report
#'    should be generated.
#' @param selectInd Numeric; Gene indexes to select for swarm optimization. If
#'  NULL the selectInd slot from the spUnsupervised object is used.
#' @param vectorize Argument to \link[pso]{psoptim}.
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
#'
NULL

#' @rdname spSwarm
#' @export

setGeneric("spSwarm", function(
  spCountsSng,
  spCountsMul,
  spUnsupervised,
  ...
){
  standardGeneric("spSwarm")
})

#' @importFrom future.apply future_lapply
#' @importFrom pso psoptim
#' @importFrom matrixStats rowSums2 rowMeans2
#' @importFrom dplyr "%>%" bind_rows
#' @importFrom purrr map
#' @rdname spSwarm
#' @export

setMethod("spSwarm", c("spCounts", "spCounts", "spUnsupervised"), function(
  spCountsSng, spCountsMul, spUnsupervised,
  maxiter = 10, swarmsize = 150, nSyntheticMultiplets = 200,
  seed = 11, norm = TRUE,
  report = FALSE, reportRate = NULL, selectInd = NULL, vectorize = FALSE,
  ...
){
    
  #put a check here to make sure all slots in the spUnsupervised object are
  #filled. This should actually be regulated by the class definition BUT you
  #should probably double check that it works as expected via unit tests.
    
  #input and input checks
  sngCPM <- getData(spCountsSng, "counts.cpm")
  mulCPM <- getData(spCountsMul, "counts.cpm")
    
  #calculate fractions
  classes <- getData(spUnsupervised, "classification")
  fractions <- rep(1.0 / length(unique(classes)), length(unique(classes)))
    
  #subset top genes for use with optimization
  #sholud also check user input selectInd
  if(is.null(selectInd)) selectInd <- getData(spUnsupervised, "selectInd")
  
  multiplets <- matrix(
    mulCPM[selectInd, ],
    ncol = ncol(mulCPM),
    dimnames = list(NULL, colnames(mulCPM))
  )
  singlets <- matrix(
    sngCPM[selectInd, ],
    ncol = ncol(sngCPM),
    dimnames = list(NULL, colnames(sngCPM))
  )
  
  #setup args for optimization
  if(report) {
    control <- list(
      maxit = maxiter, s = swarmsize, trace = 1,
      REPORT = reportRate, trace.stats = TRUE
    )
    stats <- list()
  } else {
    control <- list(maxit = maxiter, s = swarmsize, vectorize = vectorize)
    stats <- list()
  }
  
  #run optimization
  to <- if(ncol(multiplets) == 1) {to <- 1} else {to <- dim(multiplets)[2]}
  
  set.seed(seed)
  opt.out <- future_lapply(
    X = 1:to, FUN = function(i) {
      .optim.fun(
        i, fractions = fractions, multiplets = multiplets, singlets = singlets,
        classes = classes, n = nSyntheticMultiplets, control = control, ...
      )
  })
  
  #process optimization results
  result <- .processResults(
    opt.out, report, norm, stats,
    unique(classes), colnames(multiplets)
  )
  
  #create object
  new("spSwarm",
    spSwarm = result[[1]], costs = result[[2]],
    convergence = result[[3]], stats = result[[4]],
    arguments = list(maxiter = maxiter, swarmsize = swarmsize)
  )
})

.processResults <- function(result, report, norm, stats, cn, rn) {
  
  #extract swarm output
  par <- data.frame(t(sapply(result, function(j) j[[1]])))
  cost <- sapply(result, function(j) j[[2]])
  counts <- t(sapply(result, function(j) j[[3]]))
  convergence <- sapply(result, function(j) j[[4]])
  convergenceKey <- c(
    "Maximal number of function evaluations reached." = 1,
    "Maximal number of iterations reached." = 2,
    "Maximal number of restarts reached." = 3,
    "Maximal number of iterations without improvement reached." = 4
  )
  convergence <- names(convergenceKey)[match(convergence, convergenceKey)]
  if(report) stats <- lapply(result, function(x) x[[6]])

  #normalize swarm output
  if(norm) {par <- par * 1/rowSums(par)}
  colnames(par) <- sort(cn)
  rownames(par) <- rn
  
  return(list(par, cost, convergence, stats))
}

.subsetSinglets <- function(classes, singlets, n) {
  purrr::map(1:n, ~sampleSinglets(classes)) %>%
    purrr::map(., ~subsetSinglets(singlets, .x)) %>%
    purrr::map(., function(x) {rownames(x) <- 1:nrow(x); x}) %>%
    do.call("rbind", .) %>%
    .[order(as.numeric(rownames(.))), ]
}

.optim.fun <- function(
  i, fractions, multiplets, singlets, classes,
  n, control, ...
){
  oneMultiplet <- ceiling(multiplets[, i])
  singletSubset <- .subsetSinglets(classes, singlets, n)
  pso::psoptim(
    par = fractions, fn = calculateCost, oneMultiplet = oneMultiplet,
    singletSubset = singletSubset, n = n, lower = 0, upper = 1,
    control = control, ...
  )
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
#' @importFrom dplyr filter pull rowwise do ungroup full_join select bind_cols
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom tidyr separate
#' @importFrom utils combn
#' @importFrom rlang .data
#' @export

spSwarmPoisson <- function(
  spSwarm,
  edge.cutoff = 0,
  min.pval = 1,
  min.num.edges = 0,
  ...
){
  mat <- getData(spSwarm, "spSwarm")
  logic <- .fractionCutoff(mat, edge.cutoff)

  #calcluate weight
  edges <- .calculateWeight(mat, logic)

  #calculate p-value
  out <- .calculateP(edges, min.pval, min.num.edges)

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
  do(bind_cols(.data, ps(edges, .data$from, .data$to, .data$weight))) %>%
  ungroup %>%
  filter(.data$pval <= min.pval & .data$weight >= min.num.edges)
}

.calculateWeight <- function(
  mat,
  logic,
  ...
){
  homo <- paste(sort(colnames(mat)), sort(colnames(mat)), sep = "-")
  hetero <- apply(combn(colnames(mat), 2), 2, function(x) {
      paste(sort(x), collapse = "-")
  })
  totcomb <- unique(c(homo, hetero))
    
  com <- apply(logic, 1, .funx, totcomb)
  rownames(com) <- totcomb
  
  res <- apply(com, 1, sum) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "value") %>%
  setNames(c("value", "weight"))
  
  totcomb %>%
    as_tibble %>%
    separate(value, c("from", "to"), sep = "-", remove = FALSE) %>%
    full_join(res, by = "value") %>%
    select(-value) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

.funx <- function(
  row,
  totcomb
){
  pick <- names(row)[which(row)]
  if(length(pick) == 1) {
    totcomb %in% paste(pick, pick, sep = "-")
  } else {
    picked <- apply(
      combn(pick, 2), 2, function(x) {
        paste(sort(x), collapse = "-")
    })
    totcomb %in% picked
  }
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
  distFun = .dtsnCellNum,
  ...
){
  #spCounts should only include multiplets
  
  groupMeans <- getData(spUnsupervised, "groupMeans")
  selectInd <- getData(spUnsupervised, "selectInd")
  frac <- getData(spSwarm, "spSwarm")
  counts <- getData(spCounts, "counts.cpm")
  
  cellTypes <- groupMeans[selectInd, ]
  multiplets <- counts[selectInd, ]
  
  if(all(c("e", "cellNumber") %in% names(list(...)))) {
    e <- list(...)[['e']]
    cellNumber <- list(...)[['cellNumber']]
    cellNumber <- cellNumber[match(colnames(multiplets), cellNumber$sampleName), "cellNumberMedian"][[1]]
  }
  
  diff <- sapply(1:ncol(multiplets), function(x) {
    distFun(as.numeric(frac[x, ]), cellTypes, multiplets[, x], e = e, cellNumber = cellNumber[x])
  })

  colnames(diff) <- rownames(frac)
  
  if(!is.null(clusters) & !is.null(edge.cutoff)) {
    diff <- diff[,
      getMultipletsForEdge(spSwarm, edge.cutoff, clusters[1], clusters[2])
    ]
  }
  return(diff)
}

.dtsnCellNum <- function(
 fractions, cellTypes, oneMultiplet, e, cellNumber, ...
){
  if(sum(fractions) == 0) {
    return(999999999)
  }
  normFractions <- fractions / sum(fractions)
  cellTypes <- cellTypes / mean(cellTypes)
  a <- .makeSyntheticSlice(cellTypes, normFractions)
  a <- a/mean(a)
  k <- length(which(normFractions > 0))
  penalty <- .complexityPenilty(k, e, cellNumber)
  abs((oneMultiplet - a) / (a + 1)) * penalty
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
#' @importFrom utils combn

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
