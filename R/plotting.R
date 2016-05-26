
#'@include All-classes.R
NULL

#' Plot Ercc Fraction
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name plotErccFrac
#' @rdname plotErccFrac
#' @aliases plotErccFrac
#' @param counts Counts object.
#' @param ... Additional arguments to pass on
#' @return Ercc fraction plot.
#' @author Jason T. Serviss
#' @keywords plotErccFrac
#' @examples
#'
#' #use demo data
#' data(Doublet_project_data)
#'
#' #run function
#' plotErccFrac(expData)
#'
NULL

#' @rdname plotErccFrac
#' @export

setGeneric("plotErccFrac", function(x, ...
){ standardGeneric("plotErccFrac") })


#' @rdname plotErccFrac
#' @export
#' @import ggplot2

setMethod("plotErccFrac", "spCounts", function(x, ...)
{
    counts <- getData(x, "counts")
    counts.ercc <- getData(x, "counts.ercc")
    groups <- getData(x, "sampleType")
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    d <- data.frame(sampleType = groups, frac.ercc=frac.ercc)
    ggplot(d, aes(x=sampleType, y=frac.ercc))+
    geom_jitter()
})

#' Plot Markers
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name plotMarkers
#' @rdname plotMarkers
#' @aliases plotMarkers
#' @param counts Counts object.
#' @param marker1 Marker to plot on the x-axis.
#' @param marker2 Marker to plot on the y-axis.
#' @param ... Additional arguments to pass on
#' @return Ercc fraction plot.
#' @author Jason T. Serviss
#' @keywords plotMarkers
#' @examples
#'
#' #use demo data
#' data(Doublet_project_data)
#'
#' #run function
#' plotMarkers(expData, marker1="A1BG", marker2="A1CF")
#'
NULL

#' @rdname plotMarkers
#' @export
setGeneric("plotMarkers", function(x, ...
){ standardGeneric("plotMarkers") })

#' @rdname plotMarkers
#' @export
#' @import ggplot2

setMethod("plotMarkers", "spCounts", function(
    x,
    marker1,
    marker2,
    ...
){
    counts.log <- getData(x, "counts.log")
    groups <- getData(x, "sampleType")
    d <- data.frame(
        sampleType = groups,
        marker1 = counts.log[marker1, ],
        marker2 = counts.log[marker2, ]
    )
    ggplot(d, aes(x=marker1, y=marker2, colour=sampleType))+
    geom_point()
})

