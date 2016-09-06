#' @include sp.scRNAseq-package.R
NULL

#' @rdname spCounts
#' @export
.Counts <- setClass("spCounts", representation(
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
.Counts <- setClass("spUnsupervised", representation(
    counts.log="matrix",
    dist="matrix",
    tsne="matrix",
    groupMeans="matrix",
    mclust="list"
))


#' @rdname spUnsupervised
#' @export
setMethod("getData", "spUnsupervised", function(x, n=NULL)
{
    if(class(n)=="character"){
        slot(x,n)
    }
})

#' @export
.Counts <- setClass("spSwarm", representation(
    spSwarm="data.frame",
    codedSwarm="data.frame",
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