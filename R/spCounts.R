
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
#' @param sampleType A character vector indicating the column naming scheme showing that the column contains a multuplet.
#' @param ... additional arguments to pass on
#' @return The spCounts function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords spCounts
#' @examples
#'
#' #use demo data
#' data(expData)
#' counts <- getData(expData, "counts")
#' counts.ercc <- getData(expData, "counts.ercc")
#' sampleType <- '1000102901' #indicates that the column contains a multuplet
#'
#' #run function
#' cObj <- spCounts(counts, counts.ercc, '1000102901')
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
    sampleType = '1000102901',
    ...
){
    new("spCounts",
        counts=counts,
        counts.log=.norm.log.counts(counts),
        counts.ercc=counts.ercc,
        sampleType=.sampleType(sampleType, counts)
    )
})


.norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t(apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
}

.sampleType <- function(sampleType, counts) {
    dbl <- rep("Singlet", length=ncol(counts))
    dbl[grepl(sampleType, colnames(counts))] <- "Multuplet"
    return(dbl)
}

#' Filter Cells
#'
#' Select a cutoff for bad cells based on actin expression.
#'
#' A method to filter low quality cells from the analysis.
#'
#' @name filterCells
#' @rdname filterCells
#' @aliases filterCells
#' @param x Dataset of class Counts.
#' @param quantile.cut "p" argument to qnorm.
#' @param gene.name Gene to filter on.
#' @param ... additional arguments to pass on
#' @return The filterCells function returns an object of class Counts.
#' @author Jason T. Serviss
#' @keywords filterCells
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname filterCells
#' @export

setGeneric("filterCells", function(x, ...
){ standardGeneric("filterCells") })

#' @rdname filterCells
#' @export
#' @importFrom stats median qnorm
setMethod("filterCells", "spCounts",
    function(
        x,
        quantile.cut = 0.001,
        gene.name = 'ACTB',
        ...
){
    counts.log <- getData(x, "counts.log")
    counts <- getData(x, "counts")
    counts.ercc <- getData(x, "counts.ercc")
    sampleType <- getData(x, "sampleType")
    
    good.cells <- counts.log['ACTB',] > .get.cutoff.lognorm(counts.log, quantile.cut, gene.name)
    
    x@counts <- counts[ ,good.cells]
    x@counts.log <- counts.log[ ,good.cells]
    x@counts.ercc <- counts.ercc[ ,good.cells]
    x@sampleType <- sampleType[good.cells]
    
    return(x)
})

.get.cutoff.lognorm <- function(my.counts.log, quantile.cut, gene.name) {
    cl.act <- my.counts.log[gene.name,]
    cl.act.m <- median(cl.act)
    cl.act.sd <- sqrt(sum((cl.act[cl.act > cl.act.m] - cl.act.m)^2)/(sum(cl.act  > cl.act.m)-1))
    my.cut <- qnorm(p=quantile.cut, mean=cl.act.m, sd=cl.act.sd)
    my.cut
}



