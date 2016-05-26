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
