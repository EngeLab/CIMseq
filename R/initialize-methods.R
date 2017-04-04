#'@include sp.scRNAseq-package.R
NULL

#' @rdname spCounts
#' @export

setMethod("initialize","spCounts", function(
    .Object,
    ...,
    counts,
    counts.log,
    counts.cpm,
    counts.ercc
){
    callNextMethod(
        .Object,
        ...,
        counts = counts,
        counts.log = counts.log,
        counts.cpm = counts.cpm,
        counts.ercc = counts.ercc
    )
})

#' @rdname spUnsupervised
#' @export

setMethod("initialize","spUnsupervised", function(
    .Object,
    ...,
    tsne,
    tsneMeans,
    groupMeans,
    classification,
    selectInd
){
    callNextMethod(
        .Object,
        ...,
        tsne = tsne,
        tsneMeans = tsneMeans,
        groupMeans = groupMeans,
        classification = classification,
        selectInd = selectInd
    )
})

#' @rdname spSwarm
#' @export

setMethod("initialize","spSwarm", function(
    .Object,
    ...,
    spSwarm,
    costSum,
    arguments
){
    callNextMethod(
    .Object,
    ...,
    spSwarm = spSwarm,
    costSum = costSum,
    arguments = arguments
    )
})