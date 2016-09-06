#'@include sp.scRNAseq-package.R
NULL

#' @rdname spCounts
#' @export

setMethod("initialize","spCounts", function(
    .Object,
    ...,
    counts,
    counts.log,
    counts.ercc,
    sampleType
){
    callNextMethod(
        .Object,
        ...,
        counts = counts,
        counts.log = counts.log,
        counts.ercc = counts.ercc,
        sampleType = sampleType
    )
})

#' @rdname spUnsupervised
#' @export

setMethod("initialize","spUnsupervised", function(
    .Object,
    ...,
    counts.log,
    dist,
    tsne,
    groupMeans,
    mclust
){
    callNextMethod(
        .Object,
        ...,
        counts.log = counts.log,
        dist = dist,
        tsne = tsne,
        groupMeans = groupMeans,
        mclust = mclust
    )
})

#' @rdname spSwarm
#' @export

setMethod("initialize","spSwarm", function(
    .Object,
    ...,
    spSwarm,
    arguments
){
    callNextMethod(
    .Object,
    ...,
    spSwarm = spSwarm,
    codedSwarm = codedSwarm,
    arguments = arguments
    )
})