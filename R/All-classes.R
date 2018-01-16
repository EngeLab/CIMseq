#' @include sp.scRNAseq-package.R
NULL

#####################
#                   #
#     spCounts      #
#                   #
#####################

#' @rdname spCounts
#' @export
.spCounts <- setClass("spCounts", representation(
  counts = "matrix",
  counts.log = "matrix",
  counts.cpm = "matrix",
  counts.ercc = "matrix"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname spCounts
setGeneric("getData", function(
  object,
  ...
){ standardGeneric("getData") })

#' @rdname spCounts
#' @export
setMethod("getData", "spCounts", function(
  object,
  n = NULL
){
  if(class(n) == "character"){
    slot(object, n)
  }
})

#####################
#                   #
#  spUnsupervised   #
#                   #
#####################

#' @rdname spUnsupervised
#' @export
.spUnsupervised <- setClass("spUnsupervised", representation(
  tsne = "matrix",
  tsneMeans = "data.frame",
  groupMeans = "matrix",
  classification = "character",
  uncertainty = "numeric",
  selectInd = "numeric"
))

#############
#           #
# Accessors #
#           #
#############
#slot access is achieved with for ex. classification(uObj) <- newClassification

#' @rdname spUnsupervised
#' @export
setMethod("getData", "spUnsupervised", function(
  object,
  n = NULL
){
  if(class(n) == "character"){
      slot(object, n)
  }
})

#' @rdname spUnsupervised
setGeneric("tsne", function(
  object
){
  standardGeneric("tsne")
})

#' @rdname spUnsupervised
#' @export
setMethod("tsne", "spUnsupervised", function(
  object
){
  object@tsne
})

#' @rdname spUnsupervised
setGeneric("tsneMeans", function(
  object
){
  standardGeneric("tsneMeans")
})

#' @rdname spUnsupervised
#' @export
setMethod("tsneMeans", "spUnsupervised", function(
  object
){
  object@tsneMeans
})

#' @rdname spUnsupervised
setGeneric("groupMeans", function(
  object
){
  standardGeneric("groupMeans")
})

#' @rdname spUnsupervised
#' @export
setMethod("groupMeans", "spUnsupervised", function(
  object
){
  object@groupMeans
})

#' @rdname spUnsupervised
setGeneric("classification", function(
  object
){
  standardGeneric("classification")
})

#' @rdname spUnsupervised
#' @export
setMethod("classification", "spUnsupervised", function(
  object
){
  object@classification
})

#' @rdname spUnsupervised
setGeneric("uncertainty", function(
  object
){
  standardGeneric("uncertainty")
})

#' @rdname spUnsupervised
#' @export
setMethod("uncertainty", "spUnsupervised", function(
  object
){
  object@uncertainty
})

#' @rdname spUnsupervised
setGeneric("selectInd", function(
  object
){
  standardGeneric("selectInd")
})

#' @rdname spUnsupervised
#' @export
setMethod("selectInd", "spUnsupervised", function(
  object
){
  object@selectInd
})

##############
#            #
# Replacment #
#            #
##############
#https://www.bioconductor.org/help/course-materials/2013/CSAMA2013/friday/afternoon/S4-tutorial.pdf

#' @rdname spUnsupervised
setGeneric("tsne<-", function(
  object,
  value
){
  standardGeneric("tsne<-")
})

#' @rdname spUnsupervised
#' @export
setMethod("tsne<-", "spUnsupervised", function(
  object,
  value
){
  object@tsne <- value
  if (validObject(object)) return(object)
})

#' @rdname spUnsupervised
setGeneric("tsneMeans<-", function(
  object,
  value
){
  standardGeneric("tsneMeans<-")
})

#' @rdname spUnsupervised
#' @export
setMethod("tsneMeans<-", "spUnsupervised", function(
  object,
  value
){
  object@tsneMeans <- value
  if (validObject(object)) return(object)
})

#' @rdname spUnsupervised
setGeneric("groupMeans<-", function(
  object,
  value
){
  standardGeneric("groupMeans<-")
})

#' @rdname spUnsupervised
#' @export
setMethod("groupMeans<-", "spUnsupervised", function(
  object,
  value
){
  object@groupMeans <- value
  if (validObject(object)) return(object)
})

#' @rdname spUnsupervised
setGeneric("classification<-", function(
  object,
  value
){
  standardGeneric("classification<-")
})

#' @rdname spUnsupervised
#' @export
setMethod("classification<-", "spUnsupervised", function(
  object,
  value
){
  object@classification <- value
  if (validObject(object)) return(object)
})

#' @rdname spUnsupervised
setGeneric("uncertainty<-", function(
  object,
  value
){
  standardGeneric("uncertainty<-")
})

#' @rdname spUnsupervised
#' @export
setMethod("uncertainty<-", "spUnsupervised", function(
  object,
  value
){
  object@uncertainty <- value
  if (validObject(object)) return(object)
})

#' @rdname spUnsupervised
setGeneric("selectInd<-", function(
  object,
  value
){
  standardGeneric("selectInd<-")
})

#' @rdname spUnsupervised
#' @export
setMethod("selectInd<-", "spUnsupervised", function(
  object,
  value
){
  object@selectInd <- value
  if (validObject(object)) return(object)
})

#####################
#                   #
#      spSwarm      #
#                   #
#####################

#' @rdname spSwarm
#' @export
.spSwarm <- setClass("spSwarm", representation(
  spSwarm = "data.frame",
  costs = "numeric",
  convergence = "character",
  stats = "list",
  arguments = "list"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname spSwarm
#' @export
setMethod("getData", "spSwarm", function(
  object,
  n = NULL
){
  if(class(n) == "character"){
      slot(object, n)
  }
})
