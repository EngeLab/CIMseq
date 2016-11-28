#'@include sp.scRNAseq-package.R
NULL

#' @rdname spCounts
#' @export

setMethod("initialize","spCounts", function(
    .Object,
    ...,
    counts,
    counts.cpm,
    counts.log,
    counts.ercc
){
    callNextMethod(
        .Object,
        ...,
        counts = counts,
        counts.cpm = counts.cpm,
        counts.log = counts.log,
        counts.ercc = counts.ercc,
    )
})

#' @rdname spUnsupervised
#' @export

setMethod("initialize","spUnsupervised", function(
    .Object,
    ...,
    counts,
    counts.log,
    sampleType,
    tsne,
    tsneMeans,
    groupMeans,
    classification,
    selectInd
){
    callNextMethod(
        .Object,
        ...,
        counts = counts,
        counts.log = counts.log,
        sampleType = sampleType,
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
    codedSwarm,
    spUnsupervised,
    arguments
){
    callNextMethod(
    .Object,
    ...,
    spSwarm = spSwarm,
    codedSwarm = codedSwarm,
    spUnsupervised = spUnsupervised,
    arguments = arguments
    )
})