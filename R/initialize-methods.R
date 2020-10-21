#'@include CIMseq-package.R
NULL

#' @rdname CIMseqSinglets
#' @export

setMethod("initialize","CIMseqSinglets", function(
  .Object, ..., counts, counts.log, counts.cpm, 
  dim.red, classification
){
  callNextMethod(
    .Object, ..., counts = counts, counts.log = counts.log,
    counts.cpm = counts.cpm,
    dim.red = dim.red, classification = classification
  )
})

#' @rdname CIMseqMultiplets
#' @export

setMethod("initialize","CIMseqMultiplets", function(
  .Object, ..., counts, counts.log, counts.cpm, features
){
  callNextMethod(
    .Object, ..., counts = counts, counts.log = counts.log,
    counts.cpm = counts.cpm, features = features
  )
})

#' @rdname CIMseqSwarm
#' @export

setMethod("initialize","CIMseqSwarm", function(
  .Object, ..., fractions, costs, convergence,
  stats, singletIdx, arguments
){
  callNextMethod(
    .Object, ..., fractions = fractions, costs = costs, 
    convergence = convergence, stats = stats, singletIdx = singletIdx, 
    arguments = arguments
  )
})
