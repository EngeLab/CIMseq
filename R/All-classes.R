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
    tsne="matrix"
))

#' @rdname spUnsupervised
setGeneric("getData", function(x, ...
){ standardGeneric("getData") })

#' @rdname spUnsupervised
#' @export
setMethod("getData", "spUnsupervised", function(x, n=NULL)
{
    if(class(n)=="character"){
        slot(x,n)
    }
})
