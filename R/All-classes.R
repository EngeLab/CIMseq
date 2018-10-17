#' @include sp.scRNAseq-package.R
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

setGeneric("getData", function(object, ...){
  standardGeneric("getData") 
})

#' @rdname CIMseqSinglets
#' @export

setMethod("getData", "CIMseqSinglets", function(object, n = NULL){
  if(is.character(n) & .hasSlot(object, n)){
    slt <- slot(object, n)
    if(!is.function(slt)) {
      return(slt)
    } else {
      return(slt(slot(object, "counts")))
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

setGeneric("getData<-", function(object, n, value){
  standardGeneric("getData<-") 
})

#' @rdname CIMseqSinglets
#' @export

setMethod("getData<-", "CIMseqSinglets", function(object, n = NULL, value){
  if(class(n) == "character" & .hasSlot(object, n)){
    slot(object, n) <- value
    return(object)
  }
})

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

setMethod("getData", "CIMseqMultiplets", function(object, n = NULL){
  if(is.character(n) & .hasSlot(object, n)){
    slt <- slot(object, n)
    if(!is.function(slt)) {
      return(slt)
    } else {
      return(slt(slot(object, "counts")))
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

setMethod("getData<-", "CIMseqMultiplets", function(object, n = NULL, value){
  if(class(n) == "character" & .hasSlot(object, n)){
    slot(object, n) <- value
    return(object)
  }
})

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

setMethod("getData", "CIMseqSwarm", function(object, n = NULL){
  if(class(n) == "character" & .hasSlot(object, n)){
    slot(object, n)
  }
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
