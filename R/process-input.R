
#'@include All-classes.R
NULL

#' Title
#'
#' Subtitle
#'
#' Description
#'
#' @name importCounts
#' @rdname importCounts
#' @aliases importCounts
#' @param counts Counts matrix with samples as columns and genes as rows.
#' @param ... additional arguments to pass on
#' @return The importCounts function returns an object of class counts
#' @author Jason T. Serviss
#' @keywords importCounts
#' @examples
#'
#' #use demo data
#' data(Doublet_project_data)
#'
#' #run function
#' counts.log <- norm.log.counts(counts)
#'
NULL

#' @rdname importCounts
#' @export

setGeneric("importCounts", function(counts, ...
){ standardGeneric("importCounts") })

#' @rdname importCounts
#' @export
setMethod("importCounts", "matrix",
function(
    counts,
    counts.ercc,
    multID = '1000102901',
    ...
){
    new("Counts",
        counts=counts,
        counts.log=.norm.log.counts(counts),
        counts.ercc=counts.ercc,
        multID=.sampleType(multID)
    )
})


.norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t( apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
}

.sampleType <- function(multID) {
    dbl <- rep("Singlet", length=length(colnames(counts)))
    dbl[grepl(multID, colnames(counts))] <- "Doublet"
    return(dbl)
}

#' Title
#'
#' Select a cutoff for bad cells based on actin expression.
#'
#' Description
#'
#' @name filterCells
#' @rdname filterCells
#' @aliases filterCells
#' @param x Dataset of class Counts.
#' @param quantile.cut
#' @param gene.name
#' @param ... additional arguments to pass on
#' @return The filterCells function returns an object of class Counts.
#' @author Jason T. Serviss
#' @keywords filterCells
#' @examples
#'
#' #use demo data
#' data(Doublet_project_data)
#' x <- importCounts(counts=counts, counts.ercc=counts.ercc)
#'
#' #run function
#' good.cells <- filterCells(x)
#'
NULL

#' @rdname filterCells
#' @export

setGeneric("filterCells", function(x, ...
){ standardGeneric("filterCells") })

#' @rdname filterCells
#' @export
setMethod("filterCells", "Counts",
    function(
        x,
        quantile.cut = 0.001,
        gene.name = 'ACTB',
        ...
){
    counts.log <- getData(x, "counts.log")
    counts <- getData(x, "counts")
    cl.act <- counts.log[gene.name,]
    cl.act.m <- median(cl.act)
    cl.act.sd <- sqrt(sum((cl.act[cl.act > cl.act.m] - cl.act.m)^2)/(sum(cl.act  > cl.act.m)-1))
    my.cut <- qnorm(p=quantile.cut, mean=cl.act.m, sd=cl.act.sd)
    good.cells <- counts.log['ACTB',] > my.cut
    x@counts <- counts[ ,good.cells]
    x@counts.log <- counts.log[ ,good.cells]
    return(x)
})





