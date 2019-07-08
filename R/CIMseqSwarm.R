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
#' @param swarmPos matrix; Final swarm positions from \link[pso]{psoptim}
#' @param x CIMseqSwarm; A CIMseqSwarm object.
#' @param object CIMseqSwarm; A CIMseqSwarm to show.
#' @param n character; Data to extract from CIMseqSwarm object.
#' @param cacheScores logical; Use score caching optimization (experimental) 
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
#' @importFrom matrixStats rowSums2 rowMeans2
#' @importFrom dplyr "%>%" bind_rows mutate
#' @importFrom purrr map map_dbl
#' @importFrom tibble tibble as_tibble add_column
#' @importFrom tidyr unnest
#' @import hashmap
#' @rdname CIMseqSwarm
#' @export

### importFrom pso psoptim

setMethod("CIMseqSwarm", c("CIMseqSinglets", "CIMseqMultiplets"), function(
  singlets, multiplets,
  maxiter = 10, swarmsize = 150, nSyntheticMultiplets = 200, seed = 11, 
  norm = TRUE, report = FALSE, reportRate = NA, vectorize = FALSE,
  permute = FALSE, singletIdx = NULL, cacheScores=FALSE, psoControl=list(),
  startSwarm = NULL, topK=NULL, ...
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
#  fractions <- rep(1.0 / length(unique(classes)), length(unique(classes)))
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
  control <- list(maxit = maxiter, s = swarmsize, vectorize = vectorize)
  if(report) {
    control <- c(control, list(
      trace = 1, REPORT = reportRate, trace.stats = TRUE
      ))
  }
  control[(namc <- names(psoControl))] <- psoControl
    
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
          if(cacheScores) {
              .optim.fun.wcache(
                  i, fractions = fractions, multiplets = mul,
                  singletSubset = t.singletSubset, n = nSyntheticMultiplets,
                  control = control, startSwarm = startSwarm, ...
              )
          } else if(!is.null(topK)){
              if(topK > length(fractions)) {
                  stop(paste("trying to use higher number of fractions than available, topK=", topK, ", num fractions=", length(fractions)))
              }
              .optim.fun.topK(
                  i, fractions = fractions, multiplets = mul,
                  singletSubset = t.singletSubset, n = nSyntheticMultiplets,
                  control = control, startSwarm = startSwarm, maxNonNull=topK, ...
              )
          } else {
              .optim.fun(
                  i, fractions = fractions, multiplets = mul,
                  singletSubset = t.singletSubset, n = nSyntheticMultiplets,
                  control = control, startSwarm = startSwarm, ...
              )
          }
      })
  
  #process and return results
  cn <- sort(unique(classes))
  rn <- colnames(mul)
    
  topKfrac <- function(xmat, maxNonNull=2) {
      t(apply(xmat, 1, function(x) {
          x[order(x, decreasing=T)[-1:-maxNonNull]] <- 0
          x/sum(x)
      }))
  }

  fractions = .processSwarm(opt.out, cn, rn, norm)
  if(!is.null(topK)) {
      fractions <- topKfrac(fractions, topK)
  }
        
  new(
    "CIMseqSwarm",
    fractions = fractions,
    costs = setNames(map_dbl(opt.out, 2), colnames(mul)),
    convergence = setNames(.processConvergence(opt.out), colnames(mul)),
    stats = if(report) {.processStats(opt.out, cn, rn)} else {tibble()},
    swarmPos = if( psoControl[['return.swarm']] ) { map(opt.out, "swarm") } else { matrix() },
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
  n, control, startSwarm = NULL, ...
){
  oneMultiplet <- round(multiplets[, i]) #change this to round() ?
  if(is.list(startSwarm)) {
      psoptim1(
          par = fractions, fn = calculateCost, oneMultiplet = oneMultiplet,
          singletSubset = singletSubset, n = n, lower = 0, upper = 1,
          control = control, Xinit = startSwarm[[i]], ...
      )
  }
  else {
      psoptim1(
          par = fractions, fn = calculateCost, oneMultiplet = oneMultiplet,
          singletSubset = singletSubset, n = n, lower = 0, upper = 1,
          control = control, Xinit = startSwarm, ...
      )
  }
}


.optim.fun.topK <- function(
  i, fractions, multiplets, singletSubset,
  n, control, startSwarm = NULL, maxNonNull=2, ...
){
  oneMultiplet <- round(multiplets[, i]) #change this to round() ?
  if(is.list(startSwarm)) {
      psoptim1(
          par = fractions, fn = calculateCostMaxNonNull, oneMultiplet = oneMultiplet,
          singletSubset = singletSubset, n = n, lower = 0, upper = 1,
          control = control, Xinit = startSwarm[[i]], maxNonNull=maxNonNull, ...
      )
  }
  else {
      psoptim1(
          par = fractions, fn = calculateCostMaxNonNull, oneMultiplet = oneMultiplet,
          singletSubset = singletSubset, n = n, lower = 0, upper = 1,
          control = control, Xinit = startSwarm, maxNonNull=maxNonNull, ...
      )
  }
}

.optim.fun.wcache <- function(
  i, fractions, multiplets, singletSubset,
  n, control, startSwarm = NULL, ...
  )
{
  oneMultiplet <- round(multiplets[, i]) #change this to int() ?
#  my.cache <- new.env(hash=TRUE)
  my.cache <- hashmap(keys="hi", values=18.2)
  psoptim1(
    par = fractions, fn = calculateCostWrapper, oneMultiplet = oneMultiplet,
    singletSubset = singletSubset, n = n, lower = 0, upper = 1,
    control = control, cache=my.cache, Xinit=startSwarm, ...
 )
}


.optim.fun.wcacheCpp <- function(
  i, fractions, multiplets, singletSubset,
  n, control, resolution=20, startSwarm = NULL, ...
  )
{
  oneMultiplet <- round(multiplets[, i]) #change this to int() ?
#  my.cache <- new.env(hash=TRUE)
  my.cache <- createHashmap();
  psoptim1(
    par = fractions, fn = calculateCostCached, oneMultiplet = oneMultiplet,
    singletSubset = singletSubset, n = n, lower = 0, upper = 1,
    control = control, cache=my.cache, resolution = resolution, Xinit = startSwarm, ...
 )
}

calculateCostWrapper <- function(oneMultiplet, singletSubset, fractions, n, cache) {
    frac.char <- paste(as.integer(fractions*20), collapse="") # Make into 0.05 increment str.
    result <- cache$find(frac.char)
#    cat(frac.char, "\t", result, "\n")
    if(is.na(result)) {
        result <- calculateCost(oneMultiplet, singletSubset, fractions, n)
        cache$insert(frac.char, result)
#        cat(frac.char, "\n")
    }
    return(result)
}

calculateCostMaxNonNull <- function(oneMultiplet, singletSubset, fractions, n, cache, maxNonNull=2) {
    fractions[order(fractions, decreasing=T)[-1:-maxNonNull]] <- 0 # Set all non-top maxNonNull fracs to 0
    calculateCost(oneMultiplet, singletSubset, fractions, n)
}

calculateCostWrapperCpp <- function(oneMultiplet, singletSubset, fractions, n, cache, resolution) {
    calculateCostCached(oneMultiplet, singletSubset, fractions, n, cache, resolution)
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
  singlets, multiplets, swarm, binary = TRUE, maxCellsPerMultiplet=Inf
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
  cnm <- estimateCells(singlets, multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet) %>%
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
  swarm, singlets, multiplets, depleted=FALSE, ...
){
  mat <- adjustFractions(singlets, multiplets, swarm, binary = TRUE, ...)

  #calcluate weight
  edges <- .calculateWeight(mat, depleted=depleted)

  #calculate p-value
  out <- .calculateP(edges, mat, depleted=depleted)

  return(out)
}

.calculateWeight <- function(mat, depleted=FALSE) {
  from <- to <- NULL
  expand.grid(
    from = colnames(mat), to = colnames(mat),
    stringsAsFactors = FALSE
  ) %>%
    filter(from != to) %>% #doesn't calculate self edges
    mutate(weight = map2_int(from, to, function(f, t) {
      sub <- mat[, colnames(mat) %in% c(f, t)]
      length(which(rowSums2(sub) == 2))
    }))
}

.calculateP <- function(
  edges, mat, depleted=FALSE, ...
){
  from <- to <- jp <- weight <- expected.edges <- NULL
  #calculate total number of edges
  total.edges <- sum(edges[, "weight"])
  edg <- edges
  #calculate expected edges
  class.freq <- colSums2(mat) #multiplet estimated cell type frequency
  names(class.freq) <- colnames(mat)
#  cs <- colSums(mat)
  
  allProbs <- expand.grid(
    from = names(class.freq), to = names(class.freq), 
    stringsAsFactors = FALSE
  ) %>%
    filter(from != to) %>%
    mutate(edges = map2_dbl(from, to, function(f, t) {
      abs <- class.freq[names(class.freq) != f]
      rel <- abs / sum(abs)
      as.numeric(rel[t])
    }))# %>%
  
  edges <- mutate(edges, frac.edges = map2_dbl(from, to, function(f, t){
      allProbs[allProbs$from == f & allProbs$to == t, "edges"]
  }))

  
  edges$totalFrom <- sapply(1:nrow(edges), function(i) {
      sum(edges$weight[edges$from == edges$from[i]])
  })
  edges$expected.edges <- edges$frac.edges * edges$totalFrom

   #   mutate(expected = sum(edg$weight))#[edg$from == from & edg$to == to])) # * edges) # Bad
  # Previously:
#    mutate(edges = map2_dbl(from, to, function(f, t) {
#      abs <- class.freq[names(class.freq) != t]
#      rel <- abs / sum(abs)
#      as.numeric(rel[f]) * class.freq[t]
#    })) %>%
#    mutate(rel = edges / sum(edges)) %>%
#    mutate(expected = total.edges * rel)

  # Previously
#  edges <- mutate(edges, expected.edges = map2_dbl(from, to, function(f, t){
#    allProbs[allProbs$from == f & allProbs$to == t, "expected"]
#  }))
  
  #calculate p-value based on observed (weight) vs. expected (expected.edges)
#  edges <- mutate(edges, pval = ppois(
#    q = weight, lambda = expected.edges, lower.tail = FALSE
#  ))
  edges$pval <- sapply(1:nrow(edges), function(i) {
      phyper(q=edges$weight[i], m=class.freq[edges$to[i]], n=sum(class.freq)-class.freq[edges$to[i]], k=sum(edges$weight[edges$from == edges$from[i]]), lower.tail=depleted)
#      phyper(q=edges$weight[i], m=class.freq[edges$to[i]], n=sum(class.freq)-class.freq[edges$to[i]], k=class.freq[edges$from[i]], lower.tail=F)
#      phyper(q=edges$weight[i], m=sum(edges$weight[edges$to == edges$to[i]]), n=sum(edges$weight[edges$to != edges$to[i] & edges$from != edges$to[i]])/2, k=sum(edges$weight[edges$from == edges$from[i]]), lower.tail=F) 
#      phyper(q=edges$weight[i], m=sum(edges$weight[edges$to == edges$to[i]]), n=sum(edges$weight[edges$to != edges$to[i]]), k=sum(edges$weight[edges$from == edges$from[i]]), lower.tail=F)
 })
#  edges$qval <- p.adjust(edges$pval, 'fdr')
  edges$pval <- p.adjust(edges$pval, 'fdr')
#  edges$qval.hyperg <- p.adjust(edges$p.hyperg, 'fdr')
  
  #calculate score = observed / expected
  if(depleted) {
      edges <- mutate(edges, score = expected.edges / weight)
  } else {  
      edges <- mutate(edges, score = weight / expected.edges)
  }
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
  swarm, singlets, multiplets, edges, ...
){
  
  edges <- mutate_if(edges, is.factor, as.character)
  fractions <- adjustFractions(singlets, multiplets, swarm)
  
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
#' @param ... additional arguments to pass on
#' @return Edge names.
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
#' @importFrom purrr map_dfr
#' @importFrom tibble  tibble
#' @importFrom utils combn
#' @importFrom matrixStats rowSums2

setGeneric("getEdgesForMultiplet", function(
  swarm, ...
){
  standardGeneric("getEdgesForMultiplet")
})

#' @rdname getEdgesForMultiplet
#' @export

setMethod("getEdgesForMultiplet", "CIMseqSwarm", function(
  swarm, singlets, multiplets, multipletName = NULL, maxCellsPerMultiplet=Inf, depleted=FALSE
){
  s <- calculateEdgeStats(swarm, singlets, multiplets, depleted=depleted)
  frac <- adjustFractions(singlets, multiplets, swarm, binary = TRUE, maxCellsPerMultiplet=maxCellsPerMultiplet)
  if(is.null(multipletName)) multipletName <- rownames(getData(swarm, "fractions"))
  frac <- frac[multipletName, , drop = FALSE]
  
  #don't count self connections or multiplets with all 0 adjusted fractions
  rs <- rowSums2(frac)
  frac <- frac[rs > 1, , drop = FALSE]
  
  l <- length(frac)
  if(l == 0) return(tibble(sample = multipletName, from = NA, to = NA))
  
  map_dfr(1:nrow(frac), function(i) {
    p.fracs <- colnames(frac)[frac[i, ] == 1]
    cmb <- expand.grid(p.fracs, p.fracs, stringsAsFactors = FALSE)
    cmb <- cmb[cmb[, 1] != cmb[, 2], ]
    tibble(
      sample = rep(rownames(frac)[i], nrow(cmb)),
      from = cmb[, 1], to = cmb[, 2]
    )
  })
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
#' @importFrom dplyr mutate select distinct
#' @importFrom tidyr unnest

setGeneric("getCellsForMultiplet", function(
  swarm, ...
){
  standardGeneric("getCellsForMultiplet")
})

#' @rdname getCellsForMultiplet
#' @export

setMethod("getCellsForMultiplet", "CIMseqSwarm", function(
  swarm, singlets, multiplets, multipletName = NULL, ...
){
  getEdgesForMultiplet(swarm, singlets, multiplets, multipletName) %>%
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

psoptim1 <- function (par, fn, gr = NULL, ..., lower=-1, upper=1,
                     start.state = NULL,
                     Xinit = NULL,
                     control = list()) {
  fn1 <- function(par) fn(par, ...)/p.fnscale
  mrunif <- function(n,m,lower,upper) {
    return(matrix(runif(n*m,0,1),nrow=n,ncol=m)*(upper-lower)+lower)
  }
  norm <- function(x) sqrt(sum(x*x))
  rsphere.unif <- function(n,r) {
    temp <- runif(n)
    return((runif(1,min=0,max=r)/norm(temp))*temp)
  }
  svect <- function(a,b,n,k) {
    temp <- rep(a,n)
    temp[k] <- b
    return(temp)
  }
  mrsphere.unif <- function(n,r) {
    m <- length(r)
    temp <- matrix(runif(n*m),n,m)
    return(temp%*%diag(runif(m,min=0,max=r)/apply(temp,2,norm)))
  }
  npar <- length(par)
  lower <- as.double(rep(lower, ,npar))
  upper <- as.double(rep(upper, ,npar))
  con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
              abstol = -Inf, reltol = 0, REPORT = 10,
              s = NA, k = 3, p = NA, w = 1/(2*log(2)),
              c.p = .5+log(2), c.g = .5+log(2), d = NA,
              v.max = NA, rand.order = TRUE, max.restart=Inf,
              maxit.stagnate = Inf, eps.stagnate = 1e-3,
              vectorize=FALSE, hybrid = FALSE, hybrid.control = NULL,
              trace.stats = FALSE, type = "SPSO2007", return.swarm = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ## Argument error checks
  if (any(upper==Inf | lower==-Inf))
    stop("fixed bounds must be provided")

  p.type <- pmatch(con[["type"]],c("SPSO2007","SPSO2011"))-1
  if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
  
  p.trace <- con[["trace"]]>0L # provide output on progress?
  p.fnscale <- con[["fnscale"]] # scale funcion by 1/fnscale
  p.maxit <- con[["maxit"]] # maximal number of iterations
  p.maxf <- con[["maxf"]] # maximal number of function evaluations
  p.abstol <- con[["abstol"]] # absolute tolerance for convergence
  p.reltol <- con[["reltol"]] # relative minimal tolerance for restarting
  p.report <- as.integer(con[["REPORT"]]) # output every REPORT iterations
  p.s <- ifelse(is.na(con[["s"]]),ifelse(p.type==0,floor(10+2*sqrt(npar)),40),
                con[["s"]]) # swarm size
  p.p <- ifelse(is.na(con[["p"]]),1-(1-1/p.s)^con[["k"]],con[["p"]]) # average % of informants
  p.w0 <- con[["w"]] # exploitation constant
  if (length(p.w0)>1) {
    p.w1 <- p.w0[2]
    p.w0 <- p.w0[1]
  } else {
    p.w1 <- p.w0
  }
  p.c.p <- con[["c.p"]] # local exploration constant
  p.c.g <- con[["c.g"]] # global exploration constant
  p.d <- ifelse(is.na(con[["d"]]),norm(upper-lower),con[["d"]]) # domain diameter
  p.vmax <- con[["v.max"]]*p.d # maximal velocity
  p.randorder <- as.logical(con[["rand.order"]]) # process particles in random order?
  p.maxrestart <- con[["max.restart"]] # maximal number of restarts
  p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
  p.epsstagnate <- con[["eps.stagnate"]] # Used for max.stagnate
  p.vectorize <- as.logical(con[["vectorize"]]) # vectorize?
  if (is.character(con[["hybrid"]])) {
    p.hybrid <- pmatch(con[["hybrid"]],c("off","on","improved"))-1
    if (is.na(p.hybrid)) stop("hybrid should be one of \"off\", \"on\", \"improved\"")
  } else {
    p.hybrid <- as.integer(as.logical(con[["hybrid"]])) # use local BFGS search
  }
  p.hcontrol <- con[["hybrid.control"]] # control parameters for hybrid optim
  if ("fnscale" %in% names(p.hcontrol))
    p.hcontrol["fnscale"] <- p.hcontrol["fnscale"]*p.fnscale
  else
    p.hcontrol["fnscale"] <- p.fnscale
  p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?
  p.returnswarm <- as.logical(con[["return.swarm"]]) # return final swarm?
  if (p.trace) {
    message("S=",p.s,", K=",con[["k"]],", p=",signif(p.p,4),", w0=",
            signif(p.w0,4),", w1=",
            signif(p.w1,4),", c.p=",signif(p.c.p,4),
            ", c.g=",signif(p.c.g,4))
    message("v.max=",signif(con[["v.max"]],4),
            ", d=",signif(p.d,4),", vectorize=",p.vectorize,
            ", hybrid=",c("off","on","improved")[p.hybrid+1])
    if (p.trace.stats) {
      stats.trace.it <- c()
      stats.trace.error <- c()
      stats.trace.f <- NULL
      stats.trace.x <- NULL
    }
  }
  ## Initialization
  if (p.reltol!=0) p.reltol <- p.reltol*p.d
  if (p.vectorize) {
    lowerM <- matrix(lower,nrow=npar,ncol=p.s)
    upperM <- matrix(upper,nrow=npar,ncol=p.s)
  }
  # Initialize solution matrix, create random if not supplied. MARTIN
  X <- Xinit
  if(is.null(X)) {
      X <- mrunif(npar,p.s,lower,upper)
  }
  # Check validity of user-supplied X
  if(nrow(X) != npar | ncol(X) != p.s)
       stop(paste("User-supplied swarm start state is not conformant with other arguments. nrow=", nrow(X), "npar=", npar, "ncol=", ncol(X), "p.s=", p.s))
  if (!any(is.na(par)) && all(par>=lower) && all(par<=upper)) X[,1] <- par
  # Initialize direction/speed vectors
  if (p.type==0) {
    V <- (mrunif(npar,p.s,lower,upper)-X)/2
  } else { ## p.type==1
    V <- matrix(runif(npar*p.s,min=as.vector(lower-X),max=as.vector(upper-X)),npar,p.s)
    p.c.p2 <- p.c.p/2 # precompute constants
    p.c.p3 <- p.c.p/3
    p.c.g3 <- p.c.g/3
    p.c.pg3 <- p.c.p3+p.c.g3
  }
  if (!is.na(p.vmax)) { # scale to maximal velocity
    temp <- apply(V,2,norm)
    temp <- pmin.int(temp,p.vmax)/temp
    V <- V%*%diag(temp)
  }
  f.x <- apply(X,2,fn1) # first evaluations
  stats.feval <- p.s
  P <- X
  f.p <- f.x
  P.improved <- rep(FALSE,p.s)
  i.best <- which.min(f.p)
  error <- f.p[i.best]
  init.links <- TRUE
  if (p.trace && p.report==1) {
    message("It 1: fitness=",signif(error,4))
    if (p.trace.stats) {
      stats.trace.it <- c(stats.trace.it,1)
      stats.trace.error <- c(stats.trace.error,error)
      stats.trace.f <- c(stats.trace.f,list(f.x))
      stats.trace.x <- c(stats.trace.x,list(X))
    }
  }
  ## Iterations
  stats.iter <- 1
  stats.restart <- 0
  stats.stagnate <- 0
  # Check stop condition
  while (stats.iter<p.maxit && stats.feval<p.maxf && error>p.abstol &&
         stats.restart<p.maxrestart && stats.stagnate<p.maxstagnate) {
    stats.iter <- stats.iter+1
    if (p.p!=1 && init.links) {
      links <- matrix(runif(p.s*p.s,0,1)<=p.p,p.s,p.s)
      diag(links) <- TRUE
    }
    ## The swarm moves
    if (!p.vectorize) {
      if (p.randorder) {
        index <- sample(p.s)
      } else {
        index <- 1:p.s
      }
      for (i in index) {
        if (p.p==1)
          j <- i.best
        else
          j <- which(links[,i])[which.min(f.p[links[,i]])] # best informant
        temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
        V[,i] <- temp*V[,i] # exploration tendency
        if (p.type==0) {
          V[,i] <- V[,i]+runif(npar,0,p.c.p)*(P[,i]-X[,i]) # exploitation
          if (i!=j) V[,i] <- V[,i]+runif(npar,0,p.c.g)*(P[,j]-X[,i])
        } else { # SPSO 2011
          if (i!=j)
            temp <- p.c.p3*P[,i]+p.c.g3*P[,j]-p.c.pg3*X[,i] # Gi-Xi
          else
            temp <- p.c.p2*P[,i]-p.c.p2*X[,i] # Gi-Xi for local=best
          V[,i] <- V[,i]+temp+rsphere.unif(npar,norm(temp))
        }
        if (!is.na(p.vmax)) {
          temp <- norm(V[,i])
          if (temp>p.vmax) V[,i] <- (p.vmax/temp)*V[,i]
        }
        X[,i] <- X[,i]+V[,i]
        ## Check bounds
        temp <- X[,i]<lower
        if (any(temp)) {
          X[temp,i] <- lower[temp]
          V[temp,i] <- 0
        }
        temp <- X[,i]>upper
        if (any(temp)) {
          X[temp,i] <- upper[temp]
          V[temp,i] <- 0
        }
        ## Evaluate function
        if (p.hybrid==1) {
          temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                        upper=upper,control=p.hcontrol)
          V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
          X[,i] <- temp$par
          f.x[i] <- temp$value
          stats.feval <- stats.feval+as.integer(temp$counts[1])
        } else {
          f.x[i] <- fn1(X[,i])
          stats.feval <- stats.feval+1
        }
        if (f.x[i]<f.p[i]) { # improvement
          P[,i] <- X[,i]
          f.p[i] <- f.x[i]
          if (f.p[i]<f.p[i.best]) {
            i.best <- i
            if (p.hybrid==2) {
              temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                            upper=upper,control=p.hcontrol)
              V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
              X[,i] <- temp$par
              P[,i] <- temp$par
              f.x[i] <- temp$value
              f.p[i] <- temp$value
              stats.feval <- stats.feval+as.integer(temp$counts[1])
            }
          }
        }
        if (stats.feval>=p.maxf) break
      }
    } else {
      if (p.p==1)
        j <- rep(i.best,p.s)
      else # best informant
        j <- sapply(1:p.s,function(i)
                    which(links[,i])[which.min(f.p[links[,i]])]) 
      temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
      V <- temp*V # exploration tendency
      if (p.type==0) {
        V <- V+mrunif(npar,p.s,0,p.c.p)*(P-X) # exploitation
        temp <- j!=(1:p.s)
        V[,temp] <- V[,temp]+mrunif(npar,sum(temp),0,p.c.p)*(P[,j[temp]]-X[,temp])
      } else { # SPSO 2011
        temp <- j==(1:p.s)
        temp <- P%*%diag(svect(p.c.p3,p.c.p2,p.s,temp))+
          P[,j]%*%diag(svect(p.c.g3,0,p.s,temp))-
          X%*%diag(svect(p.c.pg3,p.c.p2,p.s,temp)) # G-X
        V <- V+temp+mrsphere.unif(npar,apply(temp,2,norm))
      }
      if (!is.na(p.vmax)) {
        temp <- apply(V,2,norm)
        temp <- pmin.int(temp,p.vmax)/temp
        V <- V%*%diag(temp)
      }
      X <- X+V
      ## Check bounds
      temp <- X<lowerM
      if (any(temp)) {
        X[temp] <- lowerM[temp] 
        V[temp] <- 0
      }
      temp <- X>upperM
      if (any(temp)) {
        X[temp] <- upperM[temp]
        V[temp] <- 0
      }
      ## Evaluate function
      if (p.hybrid==1) { # not really vectorizing
        for (i in 1:p.s) {
          temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                        upper=upper,control=p.hcontrol)
          V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
          X[,i] <- temp$par
          f.x[i] <- temp$value
          stats.feval <- stats.feval+as.integer(temp$counts[1])
        }
      } else {
        f.x <- apply(X,2,fn1)
        stats.feval <- stats.feval+p.s
      }
      temp <- f.x<f.p
      if (any(temp)) { # improvement
        P[,temp] <- X[,temp]
        f.p[temp] <- f.x[temp]
        i.best <- which.min(f.p)
        if (temp[i.best] && p.hybrid==2) { # overall improvement
          temp <- optim(X[,i.best],fn,gr,...,method="L-BFGS-B",lower=lower,
                        upper=upper,control=p.hcontrol)
          V[,i.best] <- V[,i.best]+temp$par-X[,i.best] # disregards any v.max imposed
          X[,i.best] <- temp$par
          P[,i.best] <- temp$par
          f.x[i.best] <- temp$value
          f.p[i.best] <- temp$value
          stats.feval <- stats.feval+as.integer(temp$counts[1])
        }
      }
      if (stats.feval>=p.maxf) break
    }
    if (p.reltol!=0) {
      d <- X-P[,i.best]
      d <- sqrt(max(colSums(d*d)))
      if (d<p.reltol) {
        X <- mrunif(npar,p.s,lower,upper)
        V <- (mrunif(npar,p.s,lower,upper)-X)/2
        if (!is.na(p.vmax)) {
          temp <- apply(V,2,norm)
          temp <- pmin.int(temp,p.vmax)/temp
          V <- V%*%diag(temp)
        }
        stats.restart <- stats.restart+1
        if (p.trace) message("It ",stats.iter,": restarting")
      }
    }
#    init.links <- f.p[i.best]==error # if no overall improvement
    init.links <- abs(f.p[i.best]-error) < p.epsstagnate # if overall improvement < eps MARTIN
    stats.stagnate <- ifelse(init.links,stats.stagnate+1,0)
    error <- f.p[i.best]
    if (p.trace && stats.iter%%p.report==0) {
      if (p.reltol!=0) 
        message("It ",stats.iter,": fitness=",signif(error,4),
                ", swarm diam.=",signif(d,4))
      else
        message("It ",stats.iter,": fitness=",signif(error,4))
      if (p.trace.stats) {
        stats.trace.it <- c(stats.trace.it,stats.iter)
        stats.trace.error <- c(stats.trace.error,error)
        stats.trace.f <- c(stats.trace.f,list(f.x))
        stats.trace.x <- c(stats.trace.x,list(X))
      }
    }
  }
  if (error<=p.abstol) {
    msg <- "Converged"
    msgcode <- 0
  } else if (stats.feval>=p.maxf) {
    msg <- "Maximal number of function evaluations reached"
    msgcode <- 1
  } else if (stats.iter>=p.maxit) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  } else if (stats.restart>=p.maxrestart) {
    msg <- "Maximal number of restarts reached"
    msgcode <- 3
  } else {
    msg <- "Maximal number of iterations without improvement reached"
    msgcode <- 4
  }
  if (p.trace) message(msg)
  o <- list(par=P[,i.best],value=f.p[i.best],
            counts=c("function"=stats.feval,"iteration"=stats.iter,
              "restarts"=stats.restart),
            convergence=msgcode,message=msg)
  if (p.trace && p.trace.stats) o <- c(o,list(stats=list(it=stats.trace.it,
                                                error=stats.trace.error,
                                                f=stats.trace.f,
                                                x=stats.trace.x)))
  if(p.returnswarm) o <- c(o,list(swarm=X))
  return(o)
}
