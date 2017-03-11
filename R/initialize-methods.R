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
    counts.ercc,
    sampleType
){
    callNextMethod(
        .Object,
        ...,
        counts = counts,
        counts.log = counts.log,
        counts.cpm = counts.cpm,
        counts.ercc = counts.ercc,
        sampleType = sampleType
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