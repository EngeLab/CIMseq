
#'@include All-classes.R
NULL

#' @rdname plotErccFraq
#' @export
setGeneric("plotErccFraq", function(x, y, ...
){ standardGeneric("plotErccFraq") })


#' @rdname plotErccFraq
#' @export
#' import ggplot2

setMethod("plotErccFraq", c("Counts", "missing"), function(x, y, ...)
{
    counts <- getData(x, "counts")
    counts.log <- getData(x, "counts.log")
    groups <- getData(x, "multID")
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    d <- data.frame(sampleType = groups, frac.ercc=frac.ercc)
    ggplot(d, aes(x=sampleType, y=frac.ercc))+
    geom_jitter()
})

#' @rdname plotMarkers
#' @export
setGeneric("plotMarkers", function(x, y, ...
){ standardGeneric("plotMarkers") })


#' @rdname plotMarkers
#' @export
#' import ggplot2

setMethod("plotMarkers", c("Counts", "missing"), function(
    x,
    y,
    marker1,
    marker2,
    ...
){
    counts.log <- getData(x, "counts.log")
    groups <- getData(x, "multID")
    d <- data.frame(
        sampleType = groups,
        marker1 = counts.log[marker1, ],
        marker2 = counts.log[marker2, ]
    )
    ggplot(d, aes(x=marker1, y=marker2, colour=sampleType))+
    geom_point()
})

