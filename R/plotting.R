
#'@include All-classes.R
NULL

#' spPlot
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name spPlot
#' @rdname spPlot
#' @aliases spPlot
#' @param counts Counts matrix with samples as columns and genes as rows.
#' @param type Can be "ercc", "markers", or .......
#' @param markers A character vector with 2 markers to plot.
#' @param ... additional arguments to pass on
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords spPlot
#' @examples
#'
#' #use demo data
#' data(Doublet_project_data)
#'
#' #run function
#' spPlot(x, type="ercc")
#'
NULL

#' @rdname spPlot
#' @export

setGeneric("spPlot", function(counts, ...
){ standardGeneric("spPlot") })

#' @rdname spPlot
#' @export
#' @import ggplot2

setMethod("spPlot", "spCounts", function(
    counts,
    type,
    markers = NULL,
    ...
){
    if( type == "ercc" ) {
        p <- .erccPlot(counts)
        p
        return(p)
    }
    if( type == "markers" ) {
        p <- .markersPlot(counts, markers)
        p
        return(p)
    }
})

.erccPlot <- function(x) {
    counts <- getData(x, "counts")
    counts.ercc <- getData(x, "counts.ercc")
    groups <- getData(x, "sampleType")
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    d <- data.frame(sampleType = groups, frac.ercc=frac.ercc)
    
    p <- ggplot(d, aes(x=sampleType, y=frac.ercc))+
    geom_jitter()
    
    return(p)
}

.markersPlot <- function(x, markers) {
    counts.log <- getData(x, "counts.log")
    groups <- getData(x, "sampleType")
    d <- data.frame(
        sampleType = groups,
        marker1 = counts.log[markers[1], ],
        marker2 = counts.log[markers[2], ]
    )
    
    p <- ggplot(d, aes(x=marker1, y=marker2, colour=sampleType))+
    geom_point()
    
    return(p)
}
