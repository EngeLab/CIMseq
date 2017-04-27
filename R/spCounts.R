
#'@include All-classes.R
NULL

#' spCounts
#'
#' Subtitle
#'
#' Imports count, sampleType, and count.ercc data to a sp.scRNAseq object.
#'
#' @name spCounts
#' @rdname spCounts
#' @aliases spCounts
#' @param counts Counts matrix with samples as columns and genes as rows.
#' @param counts.log Log2 normalized counts per million.
#' @param counts.ercc A matrix containing ercc spike-in reads and their counts.
#' @param sampleType A character indicating the column naming scheme showing that the column contains a multuplet.
#' @param object spCounts object.
#' @param n Data to extract from spCounts object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spCounts object.
#' @author Jason T. Serviss
#' @keywords spCounts
#' @examples
#'
#' #run function
#' cObj <- spCounts(testCounts, testErcc, 'm.')
#'
NULL

#' @rdname spCounts
#' @export

setGeneric("spCounts", function(counts, ...
){ standardGeneric("spCounts") })

#' @rdname spCounts
#' @export
setMethod("spCounts", "matrix",
function(
    counts,
    counts.ercc,
    ...
){
    if((dim(counts)[2]) != (dim(counts.ercc)[2])) {
        cat("ncol(counts) != ncol(counts.ercc).")
    }

    new("spCounts",
        counts=counts,
        counts.log=.norm.log.counts(counts),
        counts.cpm=.norm.counts(counts),
        counts.ercc=counts.ercc
    )
})


.norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t(apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
}

.norm.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.cpm <- t(apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
}

.sampleType <- function(sampleType, counts) {
    dbl <- rep("Singlet", length=ncol(counts))
    dbl[grepl(sampleType, colnames(counts))] <- "Multuplet"
    return(dbl)
}


