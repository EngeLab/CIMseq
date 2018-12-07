#' @include CIMseq-package.R
NULL

#####################
#                   #
#   CIMseqSinglets  #
#                   #
#####################

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

######################
#                    #
#  CIMseqMultiplets  #
#                    #
######################

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
#https://www.bioconductor.org/help/course-materials/2013/CSAMA2013/friday/afternoon/S4-tutorial.pdf

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

#####################
#                   #
#    CIMseqSwarm    #
#                   #
#####################

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

#' @rdname CIMseqSwarm
#' @export

setMethod("c", c("CIMseqSwarm"), function(x, ...){
  objs <- c(list(x), list(...))
  #you probably want to do some checks here
  new("CIMseqSwarm",
    fractions = lapply(objs, getData, "fractions") %>% do.call("rbind", .),
    costs = lapply(objs, getData, "costs") %>% do.call("c", .),
    convergence = lapply(objs, getData, "convergence") %>% do.call("c", .),
    stats = lapply(objs, getData, "stats") %>% do.call("rbind", .),
    singletIdx = lapply(objs, getData, "singletIdx") %>% do.call("c", .),
    arguments = lapply(objs, getData, "arguments") %>% do.call("rbind", .)
  )
})


#####################
#                   #
#      ggplot2      #
#                   #
#####################

#' "gg" class
#'
#' @name gg-class
#' @aliases gg
#' @family gg
#'
#' @exportClass gg
setOldClass("gg")
