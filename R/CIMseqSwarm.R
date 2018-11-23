#'@include All-classes.R
NULL

#' CIMseqSwarm
#'
#' Subtitle
#'
#' Description
#'
#' @name CIMseqSwarm
#' @rdname CIMseqSwarm
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param maxiter integer; The maximum optimization iterations.
#' @param swarmsize integer; The number of swarm particals.
#' @param nSyntheticMultiplets Numeric; Value indicating the number of synthetic
#'  multiplets to generate during deconvolution.
#' @param seed integer; Seed to set before running.
#' @param norm logical; Indicates if the sum of fractions should equal 1.
#' @param report logical; Indicates if additional reporting from the
#'   optimization should be included.
#' @param reportRate integer; If report is TRUE, the iteration interval for 
#' which a report should be generated.
#' @param vectorize logical, Argument to \link[pso]{psoptim}.
#' @param permute logical; indicates if genes should be permuted before
#'  deconvolution. For use with permutation testing.
#'  @param singletIdx list; Singlet indexes to be used to choose singlets and
#'  synthesize synthetic multiplets. Facilitates using the same synthetic set 
#'  e.g. with repeated runs or permutation. 
#' @param fractions matrix; The deconvolution results.
#' @param costs numeric; The costs after optimization.
#' @param convergence character; The convergence output from \link[pso]{psoptim}.
#' One value per multiplet.
#' @param stats tbl_df; The stats output from \link[pso]{psoptim}
#' @param arguments list; Arguments passed to the CIMseqSwarm function.
#' @param singletIdx list; Indexes indicating singlets that were subset to 
#' synthesize synthetic multiplets. Facilitates recreation of the synthetic 
#' multiplets downstream.
#' @param x CIMseqSwarm; A CIMseqSwarm object.
#' @param object CIMseqSwarm; A CIMseqSwarm to show.
#' @param n character; Data to extract from CIMseqSwarm object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on.
#' @return CIMseqSwarm output.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#'
NULL

#' @rdname CIMseqSwarm
#' @export

setGeneric("CIMseqSwarm", function(
  singlets, multiplets, ...
){
  standardGeneric("CIMseqSwarm")
})

#' @importFrom future.apply future_lapply
#' @importFrom pso psoptim
#' @importFrom matrixStats rowSums2 rowMeans2
#' @importFrom dplyr "%>%" bind_rows mutate
#' @importFrom purrr map map_dbl
#' @importFrom tibble tibble as_tibble add_column
#' @importFrom tidyr unnest
#' @rdname CIMseqSwarm
#' @export

setMethod("CIMseqSwarm", c("CIMseqSinglets", "CIMseqMultiplets"), function(
  singlets, multiplets,
  maxiter = 10, swarmsize = 150, nSyntheticMultiplets = 200, seed = 11, 
  norm = TRUE, report = FALSE, reportRate = NA, vectorize = FALSE,
  permute = FALSE, singletIdx = NULL, ...
){
    
  #put a check here to make sure all slots in the spUnsupervised object are
  #filled. This should actually be regulated by the class definition BUT you
  #should probably double check that it works as expected via unit tests.
    
  #input and input checks
  sngCPM <- getData(singlets, "counts.cpm")
  mulCPM <- getData(multiplets, "counts.cpm")
    
  #calculate fractions
  classes <- getData(singlets, "classification")
  fractions <- rep(1.0 / length(unique(classes)), length(unique(classes)))
    
  #subset top genes for use with optimization
  #sholud also check user input selectInd
  selectInd <- getData(multiplets, "features")
  
  mul <- matrix(
    mulCPM[selectInd, ],
    ncol = ncol(mulCPM),
    dimnames = list(rownames(mulCPM)[selectInd], colnames(mulCPM))
  )
  sng <- matrix(
    sngCPM[selectInd, ],
    ncol = ncol(sngCPM),
    dimnames = list(rownames(sngCPM)[selectInd], colnames(sngCPM))
  )
  
  #setup args for optimization
  if(report) {
    control <- list(
      maxit = maxiter, s = swarmsize, trace = 1,
      REPORT = reportRate, trace.stats = TRUE
    )
  } else {
    control <- list(maxit = maxiter, s = swarmsize, vectorize = vectorize)
  }
  
  #run optimization
  to <- if(ncol(mul) == 1) {to <- 1} else {to <- dim(mul)[2]}
  
  #setup singlets for synthetic multiplets synthesis
  set.seed(seed)
  if(permute) {sng <- .permuteGenes(sng)}
  
  if(is.null(singletIdx)) {
    singletIdx <- purrr::map(1:nSyntheticMultiplets, ~sampleSinglets(classes))
  }
  
  singletSubset <- appropriateSinglets(singlets, singletIdx, selectInd)
  
  #deconvolution
  opt.out <- future_lapply(
    X = 1:to, FUN = function(i) {
      .optim.fun(
        i, fractions = fractions, multiplets = mul,
        singletSubset = singletSubset, n = nSyntheticMultiplets,
        control = control, ...
      )
  })
  
  #process and return results
  cn <- sort(unique(classes))
  rn <- colnames(mul)
  
  new("CIMseqSwarm",
    fractions = .processSwarm(opt.out, cn, rn, norm),
    costs = map_dbl(opt.out, 2),
    convergence = .processConvergence(opt.out),
    stats = if(report) {.processStats(opt.out, cn, rn)} else {tibble()},
    singletIdx = map(singletIdx, as.integer),
    arguments = tibble(
      maxiter = maxiter, swarmsize = swarmsize,
      nSyntheticMultiplets = nSyntheticMultiplets, seed = seed, norm = norm,
      report = report, reportRate = reportRate, features = list(selectInd),
      vectorize = vectorize, permute = permute
    )
  )
})

.processSwarm <- function(opt.out, cn, rn, norm) {
  par <- map(opt.out, 1) %>%
    do.call("rbind", .)
  
  if(norm) {par <- par * 1 / rowSums(par)}
  colnames(par) <- sort(cn)
  rownames(par) <- rn
  par
}

.processConvergence <- function(opt.out) {
  convergence <- map_dbl(opt.out, 4)
  convergenceKey <- c(
  "Maximal number of function evaluations reached.",
  "Maximal number of iterations reached.",
  "Maximal number of restarts reached.",
  "Maximal number of iterations without improvement reached."
  )
  convergenceKey[convergence]
}

.processStats <- function(opt.out, cn, rn) {
  position <- NULL
  stats <- map(opt.out, 6)

  tibble(
    sample = rn,
    iteration = map(stats, function(x) x$it),
    error = map(stats, function(x) x$error),
    fitness = map(stats, function(x) x$f),
    position = map(stats, function(x) {
      map(x$x, function(y) t(y) * 1/colSums(y))
    })
  ) %>%
    unnest() %>%
    mutate(position = map(position, function(x) {
      x %>%
      as.data.frame() %>%
      setNames(cn) %>%
      as_tibble() %>%
      add_column(swarmMemberID = 1:nrow(.), .before = 1)
    }))
}

.optim.fun <- function(
  i, fractions, multiplets, singletSubset,
  n, control, ...
){
  oneMultiplet <- round(multiplets[, i]) #change this to round() ?
  pso::psoptim(
    par = fractions, fn = calculateCost, oneMultiplet = oneMultiplet,
    singletSubset = singletSubset, n = n, lower = 0, upper = 1,
    control = control, ...
  )
}

.permuteGenes <- function(counts){
  t(apply(counts, 1, sample))
}

.costCalculationR <- function(oneMultiplet, syntheticMultiplets) {
  dpois <- NULL
  dpois(round(oneMultiplet), lambda = syntheticMultiplets) %>%
  matrixStats::rowMeans2() %>%
  log10() %>%
  ifelse(is.infinite(.) & . < 0, -323.0052, .) %>%
  sum() %>%
  `-` (.)
}

#' appropriateSinglets
#'
#' Subtitle
#'
#' Description
#'
#' @name appropriateSinglets
#' @rdname appropriateSinglets
#' @param singlets A CIMseqSinglets object.
#' @param idx numeric; Singlet indices to subset. Generated with the 
#' \code{\link{sampleSinglets}} function. THIS IS ZERO BASED since upstream 
#' calculations are done in C++.
#' @param features numeric; Indices of selected features used for deconvolution.
#' @param ... additional arguments to pass on
#' @return Appropriated singlets.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#'
#'
NULL

#' @rdname appropriateSinglets
#' @importFrom purrr map
#' @importFrom dplyr "%>%"
#' @export

appropriateSinglets <- function(
  singlets, idx, features
){
  classes <- getData(singlets, "classification")
  sngCPM <- getData(singlets, "counts.cpm")
  singlets <- matrix(
    sngCPM[features, ],
    ncol = ncol(sngCPM),
    dimnames = list(rownames(sngCPM)[features], colnames(sngCPM))
  )
  
  sub <- idx %>%
    purrr::map(., ~subsetSinglets(singlets, .x)) %>%
    purrr::map(., function(x) {rownames(x) <- 1:nrow(x); x}) %>%
    do.call("rbind", .) %>%
    .[order(as.numeric(rownames(.))), ]
  
  rownames(sub) <- paste(
    rep(rownames(singlets), each = length(idx)), 
    1:length(idx), sep = "."
  )
  colnames(sub) <- sort(unique(classes))
  sub
}

.backTransform <- function(singletSubset, n) {
  cn <- paste(rep(colnames(singletSubset), each = n), 1:n, sep = "_")
  genes <- str_replace(rownames(singletSubset), "(.*)\\..*", "\\1")
  rn <- parse_factor(
    genes,
    levels = unique(rn)
  )
  
  out <- split(singletSubset, rn) %>%
    map(~matrix(.x, nrow = 1, dimnames = list(NULL, cn))) %>%
    do.call("rbind", .)
    
  rownames(out) <- unique(genes)
  out
}

#' spSwarmPoisson
#'
#' Subtitle
#'
#' Description
#'
#' @name spSwarmPoisson
#' @rdname spSwarmPoisson
#' @param swarm A CIMseqSwarm object.
#' @param edge.cutoff numeric; The minimum fraction to consider (?).
#' @param min.pval numeric; Minimum p-value to report.
#' @param min.num.edges numeric; Minimum number of observed edges to report a 
#' connection.
#' @param ... additional arguments to pass on
#' @return CIMseqSwarm connection weights and p-values.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#' output <- spSwarmPoisson(CIMseqSwarm_test, 1/10.5)
#'
#'
NULL

#' @rdname spSwarmPoisson
#' @importFrom stats ppois
#' @importFrom dplyr filter pull rowwise do ungroup full_join select bind_cols mutate
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom purrr map2_dbl
#' @importFrom tidyr separate
#' @importFrom utils combn
#' @importFrom rlang .data
#' @export

spSwarmPoisson <- function(
  swarm,
  edge.cutoff = 0,
  min.pval = 1,
  min.num.edges = 0,
  ...
){
  mat <- getData(swarm, "fractions")
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

.calculateWeight <- function(
  mat, logic, ...
){
  value <- NULL
  comb <- expand.grid(colnames(mat), colnames(mat), stringsAsFactors = FALSE)
  totcomb <- paste(comb$Var1, comb$Var2, sep = "-")
  
  com <- apply(logic, 1, function(row) {
    pick <- names(row)[which(row)]
    if(length(pick) == 1) {
      totcomb %in% paste(pick, pick, sep = "-")
    } else {
      picked <- apply(
        expand.grid(pick, pick, stringsAsFactors = FALSE), 1, function(x) {
          paste(x, collapse = "-")
        }
      )
      totcomb %in% picked
    }
  })
  rownames(com) <- totcomb
  
  res <- apply(com, 1, sum) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(var = "value") %>%
    setNames(c("value", "weight"))
  
  totcomb %>%
    as_tibble() %>%
    separate(value, c("from", "to"), sep = "-", remove = FALSE) %>%
    full_join(res, by = "value") %>%
    select(-value) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

.calculateP <- function(
  edges, min.pval, min.num.edges, ...
){
  classes <- unique(c(edges$from, edges$to))
  
  #calculate total number of edges
  #don't account for self connections since they are underestimated in deconvolution
  total.edges <- sum(edges[edges$from != edges$to, "weight"]) * 2
  
  #calculate expected edges
  #don't account for self connections
  ct.freq <- as.numeric(table(classes) / length(classes))
  names(ct.freq) <- names(table(classes))
  allProbs <- expand.grid(names(ct.freq), names(ct.freq))
  allProbs$jp <- ct.freq[allProbs$Var1] * ct.freq[allProbs$Var2]
  allProbs <- allProbs[allProbs$Var1 != allProbs$Var2, ]
  allProbs$jp <- allProbs$jp / sum(allProbs$jp)
  edges <- mutate(
    edges, expected.edges = map2_dbl(from, to, function(f, t) {
      if(f == t) {
        NA #self conections underestimated in the deconvolution
      } else {
        total.edges * allProbs[allProbs$Var2 == f & allProbs$Var1 == t, "jp"]
      }
  }))
  
  #calculate p-value based on observed (weight) vs. expected (expected.edges)
  edges$pval <- ppois(
    q = edges$weight, lambda = edges$expected.edges, lower.tail = FALSE
  )
  
  #calculate score = observed / expected
  edges$score <- edges$weight / edges$expected.edges
  edges
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
#' @param multiplets A CIMseqMultiplets object.
#' @param swarm A CIMseqSwarm object.
#' @param ... additional arguments to pass on
#' @return Residuals (add more description).
#' @author Jason T. Serviss
#'
NULL

#' @rdname calcResiduals
#' @export
#' @importFrom purrr map_dfc set_names

calcResiduals <- function(
  multiplets,
  swarm,
  ...
){
  gene <- residual <- NULL
  frac <- getData(swarm, "fractions")
  sm <- getData(swarm, "syntheticMultiplets")
  n <- getData(swarm, "arguments")[['nSyntheticMultiplets']]
  selectInd <- getData(multiplets, "features")
  
  mulCPM <- getData(multiplets, "counts.cpm")
  multiplets <- matrix(
    mulCPM[selectInd, ],
    ncol = ncol(mulCPM),
    dimnames = list(rownames(mulCPM)[selectInd], colnames(mulCPM))
  )
  
  map_dfc(1:ncol(multiplets), function(i) {
    as.numeric(frac[rownames(frac) == colnames(multiplets)[i], ]) %>%
    adjustAccordingToFractions(sm) %>%
    multipletSums() %>%
    vecToMat(nrow(multiplets), n) %>% #double check that this is happening as expected
    calculateCostDensity(multiplets[, i], .) %>%
    calculateLogRowMeans() %>%
    fixNegInf() %>%
    multiply_by(-1) %>%
    matrix_to_tibble(drop = TRUE)
  }) %>%
  set_names(colnames(multiplets)) %>%
  add_column(gene = rownames(multiplets), .before = 1) %>%
  gather(sample, residual, -gene)
}

#' getMultipletsForEdge
#'
#' Returns the names of the multiplets that are associated with an edge.
#'
#' Description
#'
#' @name getMultipletsForEdge
#' @rdname getMultipletsForEdge
#' @param swarm A CIMseqSwarm object.
#' @param edge.cutoff numeric; The minimum fraction to consider (?).
#' @param edges data.frame; Edges of interest. Edges are indicated by the nodes
#' they connect with one node in column one and the other node in column 2.
#' @param ... additional arguments to pass on
#' @return If the edges argument only includes one row, a vector of multiplet
#'    names is returned. If several edges are interogated a list is returned
#'    with one element per edge containing the names of the multiplets.
#' @author Jason T. Serviss
#' @examples
#'
#' e <- data.frame("A375", "HOS")
#' output <- getMultipletsForEdge(CIMseqSwarm_test, 1/10.5, e)
#'
NULL

#' @rdname getMultipletsForEdge
#' @export

setGeneric("getMultipletsForEdge", function(
  swarm,
  ...
){
    standardGeneric("getMultipletsForEdge")
})

#' @rdname getMultipletsForEdge
#' @importFrom rlang .data
#' @importFrom dplyr mutate mutate_if select rename pull
#' @export

setMethod("getMultipletsForEdge", "CIMseqSwarm", function(
  swarm, edge.cutoff, edges, ...
){
  
  edges <- mutate_if(edges, is.factor, as.character)
  
  mulForEdges <- lapply(1:nrow(edges), function(j) {
    cols <- c(pull(edges, 1)[j], pull(edges, 2)[j])
    
    if(identical(cols[1], cols[2])) {
      .self(j, cols, edges, swarm, edge.cutoff)
    } else {
      .nonSelf(j, cols, edges, swarm, edge.cutoff)
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
    mat <- getData(spSwarm, "fractions")
    logic <- mat > edge.cutoff
    rs <- rowSums(logic)
    
    frac <- getData(spSwarm, "fractions")[, unique(cols)]
    o <- which(frac > edge.cutoff & rs == 1)
    rownames(getData(spSwarm, "fractions"))[o]
}

.nonSelf <- function(j, cols, edges, spSwarm, edge.cutoff) {
    frac <- getData(spSwarm, "fractions")[, cols]
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
#' @param swarm A CIMseqSwarm object.
#' @param edge.cutoff numeric; The minimum fraction to consider (?).
#' @param multiplet character; The name of the multiplet of interest.
#' @param ... additional arguments to pass on
#' @return Edge names.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- getEdgesForMultiplet(CIMseqSwarm_test, 1/10.5, "m.NJB00204.G04")
#'
NULL

#' @rdname getEdgesForMultiplet
#' @export
#' @importFrom rlang .data
#' @importFrom purrr map_dfr
#' @importFrom tibble  tibble
#' @importFrom utils combn

setGeneric("getEdgesForMultiplet", function(
  swarm, ...
){
  standardGeneric("getEdgesForMultiplet")
})

#' @rdname getEdgesForMultiplet
#' @export

setMethod("getEdgesForMultiplet", "CIMseqSwarm", function(
  swarm, edge.cutoff, multiplet, ...
){
  s <- spSwarmPoisson(swarm, edge.cutoff = edge.cutoff)
  frac <- getData(swarm, "fractions")[multiplet, ]
  if(length(multiplet) == 1) {
    .edgeFunSingle(multiplet, edge.cutoff, frac, s)
  } else {
    map_dfr(1:length(multiplet), function(i) {
      .edgeFunSingle(multiplet[i], edge.cutoff, frac[i, ], s)
    })
  }
})

.edgeFunSingle <- function(multiplet, edge.cutoff, frac, s) {
  keep <- frac > edge.cutoff
  n <- names(frac)[keep]
  if(length(n) == 1) {
    output <- tibble(multiplet = multiplet, from = n, to = n)
    return(output)
  } else {
    c <- combn(n, 2)
    output <- tibble(multiplet = multiplet, from = c[1, ], to = c[2, ])
    return(output)
  }
}

#' calculateCosts
#'
#'
#' Description
#'
#' @name calculateCosts
#' @rdname calculateCosts
#' @aliases calculateCosts
#' @param multiplets A CIMseqMultiplets object.
#' @param swarm A CIMseqSwarm object.
#' @param fractions WILL PROBABLY BE REMOVED
#' @param ... additional arguments to pass on
#' @return Costs
#' @author Jason T. Serviss
#' @keywords calculateCosts
#' @examples
#'
#' #
#'
NULL

#' @rdname calculateCosts
#' @export

setGeneric("calculateCosts", function(
  multiplets,
  swarm,
  fractions,
  ...
){
  standardGeneric("calculateCosts")
})

#' @rdname calculateCosts
#' @export

setMethod("calculateCosts", c("CIMseqSinglets", "CIMseqSwarm", "numeric"), function(
  multiplets,
  swarm,
  fractions = NULL,
  ...
){
  if(is.null(fractions)) fractions <- getData(swarm, "fractions")
  mulCPM <- getData(multiplets, "counts.cpm")
  selectInd <- getData(swarm, "arguments")$features
  
  multiplets <- matrix(
    mulCPM[selectInd, ],
    ncol = ncol(mulCPM),
    dimnames = list(NULL, colnames(mulCPM))
  )
  
  #run optimization
  to <- if(ncol(multiplets) == 1) {to <- 1} else {to <- dim(multiplets)[2]}
  
  #setup synthetic multiplets
  singletSubset <- getData(swarm, "syntheticMultiplets")
  n <- getData(swarm, "arguments")$nSyntheticMultiplets
  
  #calculate costs
  opt.out <- future_lapply(
    X = 1:to, FUN = function(i) {
      oneMultiplet <- ceiling(multiplets[, i])
      calculateCost(oneMultiplet, singletSubset, as.numeric(fractions[i, ]), n)
  })
  names(opt.out) <- colnames(multiplets)
  opt.out
})

