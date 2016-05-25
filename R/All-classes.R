#'@include sp.scRNAseq-package.R
NULL

#' @rdname counts
#' @export
.Counts <- setClass("Counts", representation(
    counts="matrix",
    counts.log="matrix",
    counts.ercc="matrix",
    multID="character"
))


#' @rdname counts
setGeneric("getData", function(x, ...
){ standardGeneric("getData") })

#' @rdname counts
#' @export
setMethod("getData", "Counts", function(x, n=NULL)
{
    if(class(n)=="character"){
        slot(x,n)
    }
})
