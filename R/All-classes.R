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
    counts.ercc="matrix",
    sampleType="character"
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
    counts="matrix",
    counts.log="matrix",
    sampleType="character",
    tsne="matrix",
    tsneMeans="data.frame",
    groupMeans="matrix",
    classification="character",
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
setGeneric("counts<-", function(object, value
) standardGeneric("counts<-"))

#' @rdname spUnsupervised
#' @export
setMethod("counts<-", "spUnsupervised", function(object, value)
{
    object@counts <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("counts.log<-", function(object, value
) standardGeneric("counts.log<-"))

#' @rdname spUnsupervised
#' @export
setMethod("counts.log<-", "spUnsupervised", function(object, value)
{
    object@counts.log <- value
    return(object)
    
})

#' @rdname spUnsupervised
setGeneric("sampleType<-", function(object, value
) standardGeneric("sampleType<-"))

#' @rdname spUnsupervised
#' @export
setMethod("sampleType<-", "spUnsupervised", function(object, value)
{
    object@sampleType <- value
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

#####################
#                   #
#      spSwarm      #
#                   #
#####################

#' @rdname spSwarm
#' @export
.spSwarm <- setClass("spSwarm", representation(
    spSwarm="data.frame",
    codedSwarm="data.frame",
    spUnsupervised="spUnsupervised",
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