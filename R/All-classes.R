#' @include CIMseq-package.R
NULL

################################################################################
#                                                                              #
#                             CIMseqSinglets                                   #
#                                                                              #
################################################################################

#' @rdname CIMseqSinglets
#' @export

setClass("CIMseqSinglets", representation(
  counts = "matrix",
  counts.log = "function",
  counts.cpm = "function",
  counts.ercc = "matrix",
  dim.red = "matrix",
  classification = "character"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname CIMseqSinglets
#' @export

setGeneric("getData", function(x, ...){
  standardGeneric("getData") 
})

#' @rdname CIMseqSinglets
#' @export

setMethod("getData", "CIMseqSinglets", function(x, n = NULL){
  if(is.character(n) & .hasSlot(x, n)){
    slt <- slot(x, n)
    if(!is.function(slt)) {
      return(slt)
    } else {
      return(slt(slot(x, "counts")))
    }
  }
})

###############
#             #
# Replacement #
#             #
###############
#https://www.bioconductor.org/help/course-materials/2013/CSAMA2013/friday/afternoon/S4-tutorial.pdf

#' @rdname CIMseqSinglets
#' @export

setGeneric("getData<-", function(x, n, value){
  standardGeneric("getData<-") 
})

#' @rdname CIMseqSinglets
#' @export

setMethod("getData<-", "CIMseqSinglets", function(x, n = NULL, value){
  if(class(n) == "character" & .hasSlot(x, n)){
    .checkCIMseqSingletsReplacement(x, n, value)
    slot(x, n) <- value
    return(x)
  }
})

.checkCIMseqSingletsReplacement <- function(x, n, value) {
  counts <- getData(x, "counts")
  counts.ercc <- getData(x, "counts.ercc")
  classification <- getData(x, "classification")
  dim.red <- getData(x, "dim.red")
  
  if(n == "classification" & length(counts) > 0) {
    stopifnot(length(classification) == ncol(counts))
  }
  if(n == "dim.red" & length(counts) > 0) {
    stopifnot(nrow(dim.red) == ncol(counts))
  }
  if(n == "counts.ercc" & length(counts) > 0) {
    stopifnot(ncol(counts.ercc) == ncol(counts))
  }
  if(n == "counts" & length(counts.ercc) > 0) {
    stopifnot(ncol(counts.ercc) == ncol(counts))
  }
}

#################
#               #
# Concatenation #
#               #
#################

#' @rdname CIMseqSinglets
#' @export

setMethod("c", c("CIMseqSinglets"), function(x, ...){
  objs <- c(list(x), list(...))
  counts <- lapply(objs, getData, "counts")
  counts.ercc <- lapply(objs, getData, "counts.ercc")
  dim.red <- lapply(objs, getData, "dim.red")
  classification <- lapply(objs, getData, "classification")
  
  #checks
  if(!Reduce(identical, lapply(counts, rownames))) stop("Counts rownames not identical.")
  if(!Reduce(identical, lapply(counts.ercc, rownames))) stop("Counts.ercc rownames not identical.")
  if(!Reduce(identical, lapply(dim.red, colnames))) stop("dim.red colnames not identical.")
  
  new("CIMseqSinglets",
      counts = do.call("cbind", counts),
      counts.log = .norm.log.counts,
      counts.cpm = .norm.counts,
      counts.ercc = do.call("cbind", counts.ercc),
      dim.red = do.call("rbind", dim.red),
      classification = do.call("c", classification)
  )
})

################################################################################
#                                                                              #
#                             CIMseqMultiplets                                 #
#                                                                              #
################################################################################

#' @rdname CIMseqMultiplets
#' @export

setClass("CIMseqMultiplets", representation(
  counts = "matrix",
  counts.log = "function",
  counts.cpm = "function",
  counts.ercc = "matrix",
  features = "integer"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname CIMseqMultiplets
#' @export

setMethod("getData", "CIMseqMultiplets", function(x, n = NULL){
  if(is.character(n) & .hasSlot(x, n)){
    slt <- slot(x, n)
    if(!is.function(slt)) {
      return(slt)
    } else {
      return(slt(slot(x, "counts")))
    }
  }
})

###############
#             #
# Replacement #
#             #
###############

#' @rdname CIMseqMultiplets
#' @export

setMethod("getData<-", "CIMseqMultiplets", function(x, n = NULL, value){
  if(class(n) == "character" & .hasSlot(x, n)){
    .checkCIMseqMultipletsReplacement(x, n, value)
    slot(x, n) <- value
    return(x)
  }
})

.checkCIMseqMultipletsReplacement <- function(x, n, value) {
  counts <- getData(x, "counts")
  counts.ercc <- getData(x, "counts.ercc")
  
  if(n == "counts.ercc" & length(counts) > 0) {
    stopifnot(ncol(counts.ercc) == ncol(counts))
  }
  if(n == "counts" & length(counts.ercc) > 0) {
    stopifnot(ncol(counts.ercc) == ncol(counts))
  }
}

#################
#               #
# Concatenation #
#               #
#################

#' @rdname CIMseqMultiplets
#' @export

setMethod("c", c("CIMseqMultiplets"), function(x, ...){
  objs <- c(list(x), list(...))
  counts <- lapply(objs, getData, "counts")
  counts.ercc <- lapply(objs, getData, "counts.ercc")
  features <- lapply(objs, getData, "features")
  
  #checks
  if(!Reduce(identical, lapply(counts, rownames))) stop("Counts rownames not identical.")
  if(!Reduce(identical, lapply(counts.ercc, rownames))) stop("Counts.ercc rownames not identical.")
  if(length(unique(features)) != 1) warning("Features not identical, concatenating.")
  
  new("CIMseqMultiplets",
      counts = do.call("cbind", counts),
      counts.log = .norm.log.counts,
      counts.cpm = .norm.counts,
      counts.ercc = do.call("cbind", counts.ercc),
      features = unique(do.call("c", features))
  )
})

################################################################################
#                                                                              #
#                             CIMseqSwarm                                      #
#                                                                              #
################################################################################

#' @rdname CIMseqSwarm
#' @export

setClass("CIMseqSwarm", representation(
  fractions = "matrix",
  costs = "numeric",
  convergence = "character",
  stats = "tbl_df",
  singletIdx = "list",
  arguments = "tbl_df"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname CIMseqSwarm
#' @export

setMethod("getData", "CIMseqSwarm", function(x, n = NULL){
  if(class(n) == "character" & .hasSlot(x, n)){
    slot(x, n)
  }
})

#################
#               #
# Concatenation #
#               #
#################

#' @rdname CIMseqSwarm
#' @importFrom dplyr bind_rows
#' @importFrom purrr flatten
#' @export

setMethod("c", c("CIMseqSwarm"), function(x, ...){
  objs <- c(list(x), list(...))
  #you probably want to do some checks here
  cn <- lapply(objs, function(x) colnames(getData(x, "fractions")))
  if(length(unique(cn)) != 1) stop("Classes do not match. Cannot concatenate.")
  
  frac <- lapply(objs, getData, "fractions") %>% do.call("rbind", .)
  si <- lapply(objs, getData, "singletIdx")
  if(length(unique(si)) == 1) {
    si <- si[[1]]
  } else {
    warning("Non-matching singletIdx. Concatenating.")
    si <- unique(flatten(si))
  }
  new("CIMseqSwarm",
    fractions = frac,
    costs = lapply(objs, getData, "costs") %>% do.call("c", .),
    convergence = lapply(objs, getData, "convergence") %>% do.call("c", .),
    stats = lapply(objs, getData, "stats") %>% do.call("bind_rows", .),
    singletIdx = si,
    arguments = lapply(objs, getData, "arguments") %>% do.call("bind_rows", .)
  )
})

#################
#               #
#  Subsetting   #
#               #
#################

#' @rdname CIMseqSinglets
#' @export

setGeneric("filterSwarm", function(
  x, ...
){
  standardGeneric("filterSwarm") 
})

#' @rdname CIMseqSwarm
#' @importFrom dplyr filter
#' @param subset Indicates samples to be retained. Can be a character vector of 
#' sample names, a logical vector, or integer vector indicating sample indices.
#' @export

setMethod("filterSwarm", c("CIMseqSwarm"), function(x, subset){
  samples <- rownames(getData(x, "fractions"))
  s <- subset
  if(is.character(s)) s <- samples %in% s
  if(is.integer(s) | is.numeric(s)) s <- 1:length(samples) %in% s
  s <- setNames(s, samples)
  
  #NOTE
  #stats may have swarm tracking info in which case there are multiple rows
  #per sample
  new("CIMseqSwarm",
      fractions = getData(x, "fractions")[s, ],
      costs = getData(x, "costs")[s],
      convergence = getData(x, "convergence")[s],
      stats = dplyr::filter(getData(x, "stats"), sample %in% names(s[s])),
      singletIdx = getData(x, "singletIdx"),
      arguments = getData(x, "arguments")
  )
})

################################################################################
#                                                                              #
#                             ggplot2 and ggraph                               #
#                                                                              #
################################################################################

#' "gg" class
#'
#' @name gg-class
#' @aliases gg
#' @family gg
#'
#' @exportClass gg
setOldClass("gg")

#' "ggraph" class
#'
#' @name ggraph-class
#' @aliases ggraph
#' @family ggraph
#'
#' @exportClass ggraph
setOldClass("ggraph")
