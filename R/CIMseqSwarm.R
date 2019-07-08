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
#' @param swarmInit matrix; Initiation positions for the swarm.
#' @param psoControl list; Additional arguments to pso.2.0 (psoptim) function.
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
  norm = TRUE, report = FALSE, reportRate = NA,
  permute = FALSE, singletIdx = NULL, swarmInit = NULL, psoControl = list(), ...
){
    
  #put a check here to make sure all slots in the spUnsupervised object are
  #filled. This should actually be regulated by the class definition BUT you
  #should probably double check that it works as expected via unit tests.
  
  #check for same genes in singlets counts and multiplets counts
  
  
  #input and input checks
  sngCPM <- getData(singlets, "counts.cpm")
  mulCPM <- getData(multiplets, "counts.cpm")
    
  #calculate fractions
  classes <- getData(singlets, "classification")
  fractions <- rep(NA, length(unique(classes)))
    
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
  control <- list(maxit = maxiter, s = swarmsize)
  control <- c(control, psoControl)
  if(report) {
    control <- c(control, list(
      trace = 1, REPORT = reportRate, trace.stats = TRUE
    ))
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
  t.singletSubset <- t(singletSubset)
  
  #deconvolution
  opt.out <- future_lapply(
    X = 1:to, FUN = function(i) {
      .optim.fun(
        i, fractions = fractions, multiplets = mul,
        singletSubset = t.singletSubset, n = nSyntheticMultiplets,
        control = control, swarmInit = swarmInit, ...
      )
  })
  
  #process and return results
  cn <- sort(unique(classes))
  rn <- colnames(mul)
  
  new(
    "CIMseqSwarm",
    fractions = .processSwarm(opt.out, cn, rn, norm),
    costs = setNames(map_dbl(opt.out, 2), colnames(mul)),
    convergence = setNames(.processConvergence(opt.out), colnames(mul)),
    stats = if(report) {.processStats(opt.out, cn, rn)} else {tibble()},
    singletIdx = map(singletIdx, as.integer),
    arguments = tibble(
      maxiter = maxiter, swarmsize = swarmsize,
      nSyntheticMultiplets = nSyntheticMultiplets, seed = seed, norm = norm,
      report = report, reportRate = reportRate, features = list(selectInd),
      permute = permute, psoControl = list(psoControl)
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
  n, control, swarmInit, ...
){
  oneMultiplet <- round(multiplets[, i])
  pso.2.0(
    par = fractions, fn = calculateCost, oneMultiplet = oneMultiplet,
    singletSubset = singletSubset, n = n, lower = 0, upper = 1,
    control = control, swarmInit = swarmInit, ...
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
#'  If null, all genes are used.
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
  singlets, idx, features = NULL
){
  classes <- getData(singlets, "classification")
  sngCPM <- getData(singlets, "counts.cpm")
  if(is.null(features)) features <- 1:nrow(sngCPM)
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
    levels = unique(genes)
  )
  
  out <- split(singletSubset, rn) %>%
    map(~matrix(.x, nrow = 1, dimnames = list(NULL, cn))) %>%
    do.call("rbind", .)
    
  rownames(out) <- unique(genes)
  out
}

#' adjustFractions
#'
#' Subtitle
#'
#' Description
#'
#' @name adjustFractions
#' @rdname adjustFractions
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param swarm CIMseqSwarm or matrix; A CIMseqSwarm object or a matrix of 
#' fractions.
#' @param binary logical; Indicates if adjusted fractions should be returned as
#' binary values.
#' @param theoretical.max integer; See \code{\link{estimateCells}}.
#' @param ... additional arguments to pass on
#' @return Adjusted fractions matrix.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#'
#'
NULL

#' @rdname adjustFractions
#' @importFrom tibble tibble
#' @importFrom dplyr "%>%" filter full_join group_by summarize pull
#' @importFrom stats median setNames
#' @importFrom rlang .data
#' @export

adjustFractions <- function(
  singlets, multiplets, swarm, binary = TRUE, theoretical.max = Inf, ...
){
  medianCellNumber <- sampleType <- estimatedCellNumber <- NULL
  if(!is.matrix(swarm)) {
    fractions <- getData(swarm, "fractions")
  } else {
    fractions <- swarm
  }
  
  #calculate median cell number per singlet class
  cnc <- cellNumberPerClass(singlets, multiplets) %>%
    {setNames(pull(., medianCellNumber), pull(., class))}
  
  cnc <- cnc[match(colnames(fractions), names(cnc))]
  if(!identical(names(cnc), colnames(fractions))) stop("cnc name mismatch")
  
  #calculate cell number per multiplet
  cnm <- estimateCells(singlets, multiplets, theoretical.max) %>%
    filter(sampleType == "Multiplet") %>%
    {setNames(pull(., estimatedCellNumber), pull(., sample))}
  
  cnm <- cnm[match(rownames(fractions), names(cnm))]
  if(!identical(names(cnm), rownames(fractions))) stop("cnm name mismatch")
  
  #adjust fractions
  frac.renorm <- t(t(fractions) / cnc)
  adjusted <- round(frac.renorm * cnm)
  if(binary) adjusted[adjusted > 0] <- 1
  return(adjusted)
}

#' calculateEdgeStats
#'
#' Subtitle
#'
#' Description
#'
#' @name calculateEdgeStats
#' @rdname calculateEdgeStats
#' @param swarm A CIMseqSwarm object.
#' @param singlets A CIMseqSinglets object.
#' @param multiplets A CIMseqMultiplets object.
#' @param theoretical.max integer; See \code{\link{estimateCells}}.
#' @param ... additional arguments to pass on
#' @return CIMseqSwarm connection weights and p-values.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- calculateEdgeStats(
#' CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
#' )
#'
NULL

#' @rdname calculateEdgeStats
#' @importFrom stats ppois
#' @importFrom dplyr filter mutate
#' @importFrom purrr map_int map2_dbl map2_int
#' @importFrom matrixStats rowSums2 colSums2
#' @importFrom rlang .data
#' @export

calculateEdgeStats <- function(
  swarm, singlets, multiplets, theoretical.max = Inf, ...
){
  adj <- adjustFractions(
    singlets, multiplets, swarm, binary = TRUE, 
    theoretical.max = theoretical.max
  )

  #calcluate weight
  edges <- .calculateWeight(adj)

  #calculate p-value
  out <- .calculateP(edges, adj)

  return(out)
}

.calculateWeight <- function(adj, ...) {
  from <- to <- NULL
  
  .f1 <- function(f, t, ...) {
    length(which(rowSums2(adj[, colnames(adj) %in% c(f, t)]) == 2))
  }
  
  expand.grid(
    from = colnames(adj), to = colnames(adj),
    stringsAsFactors = FALSE
  ) %>%
    filter(from != to) %>% #doesn't calculate self edges
    # https://github.com/r-lib/covr/issues/377
    # mutate(weight = map2_int(from, to, function(f, t) {
    #   sub <- mat[, colnames(mat) %in% c(f, t)]
    #   length(which(rowSums2(sub) == 2))
    # }))
    mutate(weight = map2_int(from, to, function(f, t) .f1(f, t)))
}

.calculateP <- function(
  edges, mat, ...
){
  from <- to <- jp <- weight <- expected.edges <- NULL
  
  #calculate expected edges
  class.freq <- colSums2(mat) #multiplet estimated cell type frequency
  names(class.freq) <- colnames(mat)
  
  # https://github.com/r-lib/covr/issues/377
  .f1 <- function(f, d) {
    freq <- class.freq[names(class.freq) != f]
    rel <- freq / sum(freq)
    rel[pull(d, to)]
  }
  
  edges <- edges %>%
    nest(-from) %>%
    mutate(to.freq = map2(from, data, ~.f1(.x, .y))) %>%
    mutate(expected.edges = map2(to.freq, data, ~sum(pull(.y, weight)) * .x)) %>%
    unnest() %>%
    select(from, to, weight, to.freq, expected.edges)
  
  #calculate score = observed / expected
  edges <- mutate(edges, score = weight / expected.edges)
  
  #calculate p-value
  edges %>%
    mutate(pval = map_dbl(1:nrow(.),
      ~phyper(
        q = edges$weight[.x], 
        m = class.freq[edges$to[.x]], 
        n = sum(class.freq) - class.freq[edges$to[.x]], 
        k = sum(edges$weight[edges$from == edges$from[.x]]), 
        lower.tail = FALSE
      )
    )) %>% 
    mutate(pval = p.adjust(pval, 'fdr'))
}

.calculateP.poisson <- function(
  edges, mat, ...
){
  from <- to <- jp <- weight <- expected.edges <- NULL
  #calculate total number of edges
  total.edges <- sum(edges[, "weight"])
  
  #calculate expected edges
  class.freq <- colSums2(mat) #multiplet estimated cell type frequency
  names(class.freq) <- colnames(mat)
  
  .f1 <- function(f, t) {
    abs <- class.freq[names(class.freq) != t]
    rel <- abs / sum(abs)
    as.numeric(rel[f]) * class.freq[t]
  }
  
  allProbs <- expand.grid(
    from = names(class.freq), to = names(class.freq), 
    stringsAsFactors = FALSE
  ) %>%
    filter(from != to) %>%
    # https://github.com/r-lib/covr/issues/377
    # mutate(edges = map2_dbl(from, to, function(f, t) {
    #   abs <- class.freq[names(class.freq) != t]
    #   rel <- abs / sum(abs)
    #   as.numeric(rel[f]) * class.freq[t]
    # })) %>%
    mutate(edges = map2_dbl(from, to, ~.f1(.x, .y))) %>%
    mutate(rel = edges / sum(edges)) %>%
    mutate(expected = total.edges * rel)
  
  # https://github.com/r-lib/covr/issues/377
  # edges <- mutate(edges, expected.edges = map2_dbl(from, to, function(f, t){
  #   allProbs[allProbs$from == f & allProbs$to == t, "expected"]
  # }))
  
  .f2 <- function(f, t) {
    allProbs[allProbs$from == f & allProbs$to == t, "expected"]
  }
  edges <- mutate(edges, expected.edges = map2_dbl(from, to, ~.f2(.x, .y)))
  
  #calculate p-value based on observed (weight) vs. expected (expected.edges)
  edges <- mutate(edges, pval = ppois(
    q = weight, lambda = expected.edges, lower.tail = FALSE
  ))
  
  #calculate score = observed / expected
  edges <- mutate(edges, score = weight / expected.edges)
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
#' @param singlets A CIMseqSinglets object.
#' @param multiplets A CIMseqMultiplets object.
#' @param swarm A CIMseqSwarm object.
#' @param include character;  If residuals should only be calculated for a 
#' subset of the multiplets, include their names here. Default is to calculate 
#' for all multiplets.
#' @param fractions matrix; A matrix of fractions. By default the fractions in 
#' the CIMseqSwarm object are used.
#' @param ... additional arguments to pass on
#' @return Residuals (add more description).
#' @author Jason T. Serviss
#'
NULL

#' @rdname calcResiduals
#' @export
#' @importFrom purrr reduce set_names
#' @importFrom future.apply future_lapply
#' @importFrom dplyr bind_cols
#' @importFrom tibble add_column
#' @importFrom tidyr gather

calcResiduals <- function(
  singlets, multiplets, swarm, include = NULL, fractions = NULL, ...
){
  residual <- gene <- NULL
  if(is.null(fractions)) {
    frac <- getData(swarm, "fractions") 
  } else {
    frac <- fractions
  }
  selectInd <- getData(multiplets, "features")
  n <- getData(swarm, "arguments")[['nSyntheticMultiplets']]
  idx <- getData(swarm, "singletIdx")
  sm <- appropriateSinglets(singlets, idx, selectInd)
  
  mulCPM <- getData(multiplets, "counts.cpm")
  if(!is.null(include) & length(include) > 1) mulCPM <- mulCPM[, include]
  if(!is.null(include) & length(include) == 1) {
    mulCPM <- matrix(
      mulCPM[, include], 
      nrow = nrow(mulCPM), 
      dimnames = list(rownames(mulCPM), include))
  }
  
  multiplets <- matrix(
    mulCPM[selectInd, ],
    ncol = ncol(mulCPM),
    dimnames = list(rownames(mulCPM)[selectInd], colnames(mulCPM))
  )
  
  future_lapply(
    X = 1:ncol(multiplets), FUN = function(i) {
      as.numeric(frac[rownames(frac) == colnames(multiplets)[i], ]) %>%
        adjustAccordingToFractions(., sm) %>%
        multipletSums() %>%
        vecToMat(nrow(multiplets), n) %>% #double check that this is happening as expected
        calculateCostDensity(round(multiplets[, i]), .) %>%
        calculateLogRowMeans() %>%
        fixNegInf() %>%
        multiply_by(-1) %>%
        matrix_to_tibble(drop = TRUE)
  }) %>%
    reduce(., bind_cols) %>%
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
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param edges data.frame; Edges of interest. Edges are indicated by the nodes
#' they connect with one node in column one and the other node in column 2.
#' @param theoretical.max integer; See \code{\link{estimateCells}}.
#' @param ... additional arguments to pass on
#' @return If the edges argument only includes one row, a vector of multiplet
#'    names is returned. If several edges are interogated a list is returned
#'    with one element per edge containing the names of the multiplets.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- getMultipletsForEdge(
#' CIMseqSwarm_test, 
#' CIMseqSinglets_test, 
#' CIMseqMultiplets_test, 
#' data.frame("A375", "HOS")
#' )
#'
NULL

#' @rdname getMultipletsForEdge
#' @export

setGeneric("getMultipletsForEdge", function(
  swarm, ...
){
  standardGeneric("getMultipletsForEdge")
})

#' @rdname getMultipletsForEdge
#' @importFrom rlang .data
#' @importFrom dplyr mutate_if
#' @importFrom purrr map_dfr
#' @importFrom matrixStats rowSums2
#' @importFrom tibble tibble
#' @export

setMethod("getMultipletsForEdge", "CIMseqSwarm", function(
  swarm, singlets, multiplets, edges, theoretical.max = Inf, ...
){
  
  edges <- mutate_if(edges, is.factor, as.character)
  fractions <- adjustFractions(
    singlets, multiplets, swarm, theoretical.max = theoretical.max
  )
  
  map_dfr(1:nrow(edges), function(i) {
    e <- as.character(edges[i, ])
    sub <- fractions[, e]
    rs <- matrixStats::rowSums2(sub)
    multiplets <- rownames(sub)[rs == 2]
    tibble(
      sample = multiplets,
      from = rep(e[1], length(multiplets)),
      to = rep(e[2], length(multiplets))
    )
  })
})

#' getEdgesForMultiplet
#'
#' Returns the names of the edges detected in a multiplet.
#'
#' Description
#'
#' @name getEdgesForMultiplet
#' @rdname getEdgesForMultiplet
#' @aliases getEdgesForMultiplet
#' @param swarm A CIMseqSwarm object.
#' @param singlets A CIMseqSinglets object.
#' @param multiplets A CIMseqMultiplets object.
#' @param multipletName character; The name of the multiplet of interest.
#' @param theoretical.max integer; See \code{\link{estimateCells}}.
#' @param drop logical; Remove self connections?
#' @param ... additional arguments to pass on
#' @return Edge names. Note that multiplets that contain no connections are not 
#'  included in the output and neither are self connections.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- getEdgesForMultiplet(
#' CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test,
#' "m.NJB00204.G04"
#' )
#'
NULL

#' @rdname getEdgesForMultiplet
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr filter everything mutate select
#' @importFrom purrr map2
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr unnest

setGeneric("getEdgesForMultiplet", function(
  swarm, ...
){
  standardGeneric("getEdgesForMultiplet")
})

#' @rdname getEdgesForMultiplet
#' @export

setMethod("getEdgesForMultiplet", "CIMseqSwarm", function(
  swarm, singlets, multiplets, multipletName = NULL, theoretical.max = Inf, 
  drop = TRUE, ...
){
  from <- to <- NULL
  if(is.null(multipletName)) multipletName <- rownames(getData(swarm, "fractions"))
  
  s <- calculateEdgeStats(
    swarm, singlets, multiplets, theoretical.max = theoretical.max
  )
  mat <- adjustFractions(
    singlets, multiplets, swarm, binary = TRUE, 
    theoretical.max = theoretical.max
  )
  
  edges <- expand.grid(
    from = colnames(mat), to = colnames(mat),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() 
  
  if(drop) edges <- filter(edges, from != to)
  if(!drop) rs <- rowSums2(mat)
  
  .f1 <- function(f, t) {
    if(f == t) {
      rownames(mat)[mat[, colnames(mat) == f] == 1 & rs == 1]
    } else {
      sub <- mat[, colnames(mat) %in% c(f, t)]
      rownames(mat)[which(rowSums2(sub) == 2)]
    }
  }
  
  data <- edges %>%
    # https://github.com/r-lib/covr/issues/377
    # mutate(sample = map2(from, to, function(f, t) {
    #   if(f == t) {
    #     rownames(mat)[mat[, colnames(mat) == f] == 1 & rs == 1]
    #   } else {
    #     sub <- mat[, colnames(mat) %in% c(f, t)]
    #     rownames(mat)[which(rowSums2(sub) == 2)]
    #   }
    # })) %>%
    mutate(sample = map2(from, to, ~.f1(.x, .y))) %>%
    unnest() %>%
    filter(sample %in% multipletName) %>%
    select(sample, everything())
  
  if(nrow(data) == 0) return(tibble(sample = multipletName, from = NA, to = NA))

  return(data)
})

#' getCellsForMultiplet
#'
#' Returns the names of the cells detected in a multiplet.
#'
#' Description
#'
#' @name getCellsForMultiplet
#' @rdname getCellsForMultiplet
#' @aliases getCellsForMultiplet
#' @param swarm A CIMseqSwarm object.
#' @param singlets A CIMseqSinglets object.
#' @param multiplets A CIMseqMultiplets object.
#' @param multipletName character; The name of the multiplet of interest.
#' @param drop logical; Remove self connections?
#' @param ... additional arguments to pass on
#' @return Edge names.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- getCellsForMultiplet(
#' CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test,
#' "m.NJB00204.G04"
#' )
#'
NULL

#' @rdname getCellsForMultiplet
#' @export
#' @importFrom rlang .data
#' @importFrom purrr map2
#' @importFrom dplyr mutate select distinct pull bind_rows
#' @importFrom tidyr unnest
#' @importFrom tibble tibble

setGeneric("getCellsForMultiplet", function(
  swarm, ...
){
  standardGeneric("getCellsForMultiplet")
})

#' @rdname getCellsForMultiplet
#' @export

setMethod("getCellsForMultiplet", "CIMseqSwarm", function(
  swarm, singlets, multiplets, multipletName = NULL, drop = TRUE, ...
){
  if(is.null(multipletName)) multipletName <- colnames(getData(multiplets, "counts"))
  
  getEdgesForMultiplet(
    swarm, singlets, multiplets, multipletName, drop = drop
  ) %>%
    mutate(cells = map2(.data$from, .data$to, ~c(.x, .y))) %>%
    select(-.data$from, -.data$to) %>%
    unnest() %>%
    distinct()
})

#' calculateCosts
#'
#'
#' Description
#'
#' @name calculateCosts
#' @rdname calculateCosts
#' @aliases calculateCosts
#' @param singlets A CIMseqSinglets object.
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
  singlets,
  multiplets,
  swarm,
  ...
){
  standardGeneric("calculateCosts")
})

#' @rdname calculateCosts
#' @export

setMethod(
  "calculateCosts", c("CIMseqSinglets", "CIMseqMultiplets", "CIMseqSwarm"), 
  function(
    singlets, multiplets, swarm, fractions = NULL, ...
){
  if(is.null(fractions)) fractions <- getData(swarm, "fractions")
  if(is.null(dim(fractions))) fractions <- matrix(fractions, ncol = length(fractions))
  mulCPM <- getData(multiplets, "counts.cpm")
  selectInd <- getData(swarm, "arguments")$features[[1]]
  
  multiplets <- matrix(
    mulCPM[selectInd, ],
    ncol = ncol(mulCPM),
    dimnames = list(NULL, colnames(mulCPM))
  )
  
  #run optimization
  to <- if(ncol(multiplets) == 1) {to <- 1} else {to <- dim(multiplets)[2]}
  
  #setup synthetic multiplets
  sngIdx <- getData(swarm, "singletIdx")
  sngSubset <- appropriateSinglets(singlets, sngIdx, selectInd)
  nSynthMul <- getData(swarm, "arguments")$nSyntheticMultiplets[[1]]
  
  #calculate costs
  opt.out <- future_lapply(
    X = 1:to, FUN = function(i) {
      oneMultiplet <- ceiling(multiplets[, i])
      calculateCost(oneMultiplet, sngSubset, as.numeric(fractions[i, ]), nSynthMul)
  })
  names(opt.out) <- colnames(multiplets)
  opt.out
})

