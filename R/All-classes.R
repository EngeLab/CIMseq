#' @include sp.scRNAseq-package.R
NULL

#' @rdname spCounts
#' @export
.spCounts <- setClass("spCounts", representation(
    counts="matrix",
    counts.log="matrix",
    counts.ercc="matrix",
    sampleType="character"
))

#' @rdname spCounts
setGeneric("getData", function(x, ...
){ standardGeneric("getData") })

#' @rdname spCounts
#' @export
setMethod("getData", "spCounts", function(x, n=NULL)
{
    if(class(n)=="character"){
        slot(x,n)
    }
})


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

#' @rdname spUnsupervised
#' @export
setMethod("getData", "spUnsupervised", function(x, n=NULL)
{
    if(class(n)=="character"){
        slot(x,n)
    }
})

#' @rdname spSwarm
#' @export
.spSwarm <- setClass("spSwarm", representation(
    spSwarm="data.frame",
    codedSwarm="data.frame",
    spUnsupervised="spUnsupervised",
    arguments="list"
))


#' @rdname spSwarm
#' @export
setMethod("getData", "spSwarm", function(x, n=NULL)
{
    if(class(n)=="character"){
        slot(x,n)
    }
})