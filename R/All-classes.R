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
    counts.cpm="matrix",
    counts.log="matrix",
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
    unsupervisedC="matrix",
    groupMeans="matrix",
    classification="character",
    selectIdx="numeric"
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
setGeneric("unsupervisedC", function(object
) standardGeneric("unsupervisedC"))

#' @rdname spUnsupervised
#' @export
setMethod("unsupervisedC", "spUnsupervised", function(object)
{
    object@unsupervisedC
})

#' @rdname spUnsupervised
setGeneric("groupMeans", function(object
) standardGeneric("groupMeans"))

#' @rdname spUnsupervised
#' @export
setMethod("groupMeans", "spUnsupervised", function(object)
{
    object@groupMeans
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

#' @rdname spUnsupervised
setGeneric("selectIdx", function(object
) standardGeneric("selectIdx"))

#' @rdname spUnsupervised
#' @export
setMethod("selectIdx", "spUnsupervised", function(object)
{
    object@selectIdx
})

##############
#            #
# Replacment #
#            #
##############

#' @rdname spUnsupervised
setGeneric("unsupervisedC<-", function(object, value
) standardGeneric("unsupervisedC<-"))

#' @rdname spUnsupervised
#' @export
setMethod("unsupervisedC<-", "spUnsupervised", function(object, value)
{
    object@unsupervisedC <- value
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
setGeneric("selectIdx<-", function(object, value
) standardGeneric("selectIdx<-"))

#' @rdname spUnsupervised
#' @export
setMethod("selectIdx<-", "spUnsupervised", function(object, value)
{
    object@selectIdx <- value
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