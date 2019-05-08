#'@include CIMseq-package.R
NULL

#' @rdname CIMseqSinglets
#' @export

setMethod("initialize","CIMseqSinglets", function(
  .Object, ..., counts, counts.log, counts.cpm, counts.ercc, 
  dim.red, classification
){
  callNextMethod(
    .Object, ..., counts = counts, counts.log = counts.log,
    counts.cpm = counts.cpm, counts.ercc = counts.ercc,
    dim.red = dim.red, classification = classification
  )
})

#' @rdname CIMseqMultiplets
#' @export

setMethod("initialize","CIMseqMultiplets", function(
  .Object, ..., counts, counts.log, counts.cpm, counts.ercc, features
){
  callNextMethod(
    .Object, ..., counts = counts, counts.log = counts.log,
    counts.cpm = counts.cpm, counts.ercc = counts.ercc, features = features
  )
})

#' @rdname CIMseqSwarm
#' @export

setMethod("initialize","CIMseqSwarm", function(
  .Object, ..., fractions, costs, convergence,
  stats, singletIdx, swarmPositions, arguments
){
  callNextMethod(
    .Object, ..., fractions = fractions, costs = costs, 
    convergence = convergence, stats = stats, singletIdx = singletIdx, 
    swarmPositions = swarmPositions, arguments = arguments
  )
})
