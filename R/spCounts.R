
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
#' @param counts.ercc A matrix containing ercc spike-in reads and their counts.
#' @param counts.log Log2 normalized counts per million.
#' @param object spCounts object.
#' @param x Default plot param, an spCounts object containing singlets.
#' @param y Default plot param, an spCounts object containing multuplets.
#' @param type Character; The type of plot desired. Currently \emph{markers} or \emph{ercc}.
#' @param markers Markers/genes to plot. Limited to two.
#' @param n Data to extract from spCounts object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spCounts object.
#' @author Jason T. Serviss
#' @keywords spCounts
#' @examples
#'
#' #run function
#' cObj <- spCounts(testCounts, testErcc)
#'
NULL

#' @rdname spCounts
#' @export

setGeneric("spCounts", function(counts, ...
){ standardGeneric("spCounts") })

#' @rdname spCounts
#' @export
setMethod("spCounts", "matrix", function(
    counts,
    counts.ercc,
    ...
){
    .inputCheckCounts(counts, counts.ercc)
    new("spCounts",
        counts=counts,
        counts.log=.norm.log.counts(counts),
        counts.cpm=.norm.counts(counts),
        counts.ercc=counts.ercc
    )
})

.inputCheckCounts <- function(counts, counts.ercc) {
    if((dim(counts)[2]) != (dim(counts.ercc)[2])) {
        message("ncol(counts) != ncol(counts.ercc).")
    }
    if(any(is.na(c(counts, counts.ercc)))) {
        message("is.na(c(counts, counts.ercc) returned TRUE")
    }
}

.norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t(apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
    return(counts.log)
}

.norm.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.cpm <- t(apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
}



#' estimateCells
#'
#' Subtitle
#'
#' Uses ERCC data to calculate the fraction of ERCC reads in the samples. In
#' addition, this function utilizes ERCC data to estimate the cell number
#' in each sample.
#'
#' @name estimateCells
#' @rdname estimateCells
#' @aliases estimateCells
#' @param spCountsSng A spCounts object with singlets.
#' @param spCountsMul A spCounts object with multiplets.
#' @param ... additional arguments to pass on
#' @return A data frame including the fraction of ercc reads and cell counts for
#'    each sample.
#' @author Jason T. Serviss
#' @keywords spCounts
#' @examples
#'
#' #run function
#'
#'
NULL

#' @export

setGeneric("estimateCells", function(
    spCountsSng,
    spCountsMul,
    ...
){
    standardGeneric("estimateCells")
})

#' @rdname estimateCells
#' @export
setMethod("estimateCells", "spCounts", function(
    spCountsSng,
    spCountsMul,
    ...
){
    sampleType <- c(
        rep(
            "Singlet",
            ncol(getData(spCountsSng, "counts"))
        ),
        rep(
            "Multiplet",
            ncol(getData(spCountsMul, "counts"))
        )
    )
    
    counts <- cbind(
        getData(spCountsSng, "counts"),
        getData(spCountsMul, "counts")
    )
    counts.ercc <- cbind(
        getData(spCountsSng, "counts.ercc"),
        getData(spCountsMul, "counts.ercc")
    )
    
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    cellNumberMedian <- median(frac.ercc[sampleType == "Singlet"]) / frac.ercc
    cellNumberMin <- quantile(frac.ercc[sampleType == "Singlet"])[2] / frac.ercc
    cellNumberMax <- quantile(frac.ercc[sampleType == "Singlet"])[4] / frac.ercc

    d <- data.frame(
        sampleType = factor(sampleType, levels=c("Singlet", "Multiplet")),
        frac.ercc = frac.ercc,
        cellNumberMin = cellNumberMin,
        cellNumberMedian = cellNumberMedian,
        cellNumberMax = cellNumberMax
    )
    return(d)
    
})
