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
        counts.ercc = counts.ercc
    )
})

#' @rdname spUnsupervised
#' @export

setMethod("initialize","spUnsupervised", function(
    .Object,
    ...,
    unsupervisedC,
    groupMeans,
    classification,
    selectIdx
){
    callNextMethod(
        .Object,
        ...,
        unsupervisedC = unsupervisedC,
        groupMeans = groupMeans,
        classification = classification,
        selectIdx = selectIdx
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