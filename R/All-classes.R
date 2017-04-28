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
    counts="matrix",
    counts.log="matrix",
    counts.cpm="matrix",
    counts.ercc="matrix"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname spCounts
setGeneric("getData", function(object, ...
){ standardGeneric("getData") })

#' @rdname spCounts
#' @export
setMethod("getData", "spCounts", function(object, n=NULL)
{
    if(class(n)=="character"){
        slot(object,n)
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
    tsne="matrix",
    tsneMeans="data.frame",
    groupMeans="matrix",
    classification="character",
    uncertainty="numeric",
    selectInd="numeric"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname spUnsupervised
#' @export
setMethod("getData", "spUnsupervised", function(object, n=NULL)
{
    if(class(n)=="character"){
        slot(object,n)
    }
})

#' @rdname spUnsupervised
setGeneric("counts", function(object
) standardGeneric("counts"))

#' @rdname spUnsupervised
#' @export
setMethod("counts", "spUnsupervised", function(object)
{
    object@counts
})

#' @rdname spUnsupervised
setGeneric("counts.log", function(object
) standardGeneric("counts.log"))

#' @rdname spUnsupervised
#' @export
setMethod("counts.log", "spUnsupervised", function(object)
{
    object@counts.log
})

#' @rdname spUnsupervised
setGeneric("sampleType", function(object
) standardGeneric("sampleType"))

#' @rdname spUnsupervised
#' @export
setMethod("sampleType", "spUnsupervised", function(object)
{
    object@sampleType
})

#' @rdname spUnsupervised
setGeneric("classification", function(object
) standardGeneric("classification"))

#' @rdname spUnsupervised
#' @export
setMethod("classification", "spUnsupervised", function(object)
{
    object@classification
})

##############
#            #
# Replacment #
#            #
##############

#' @rdname spUnsupervised
setGeneric("tsne<-", function(object, value
) standardGeneric("tsne<-"))

#' @rdname spUnsupervised
#' @export
setMethod("tsne<-", "spUnsupervised", function(object, value)
{
    object@tsne <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("tsneMeans<-", function(object, value
) standardGeneric("tsneMeans<-"))

#' @rdname spUnsupervised
#' @export
setMethod("tsneMeans<-", "spUnsupervised", function(object, value)
{
    object@tsneMeans <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("groupMeans<-", function(object, value
) standardGeneric("groupMeans<-"))

#' @rdname spUnsupervised
#' @export
setMethod("groupMeans<-", "spUnsupervised", function(object, value)
{
    object@groupMeans <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("classification<-", function(object, value
) standardGeneric("classification<-"))

#' @rdname spUnsupervised
#' @export
setMethod("classification<-", "spUnsupervised", function(object, value)
{
    object@classification <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("uncertainty<-", function(object, value
) standardGeneric("uncertainty<-"))

#' @rdname spUnsupervised
#' @export
setMethod("uncertainty<-", "spUnsupervised", function(object, value)
{
    object@uncertainty <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("selectInd<-", function(object, value
) standardGeneric("selectInd<-"))

#' @rdname spUnsupervised
#' @export
setMethod("selectInd<-", "spUnsupervised", function(object, value)
{
    object@selectInd <- value
    return(object)
    
})

#####################
#                   #
#      spSwarm      #
#                   #
#####################

#' @rdname spSwarm
#' @export
.spSwarm <- setClass("spSwarm", representation(
    spSwarm="data.frame",
    costSum="numeric",
    arguments="list"
))

#############
#           #
# Accessors #
#           #
#############

#' @rdname spSwarm
#' @export
setMethod("getData", "spSwarm", function(object, n=NULL)
{
    if(class(n)=="character"){
        slot(object,n)
    }
})