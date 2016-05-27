
#'@include All-classes.R
NULL

#' Run Tsne
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name run.tsne
#' @rdname run.tsne
#' @aliases run.tsne
#' @param counts Counts object.
#' @param ... Additional arguments to pass on
#' @return Ercc fraction plot.
#' @author Jason T. Serviss
#' @keywords run.tsne
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname run.tsne
#' @export

setGeneric("run.tsne", function(x, ...
){ standardGeneric("run.tsne") })


#' @rdname run.tsne
#' @export
#' @importFrom tsne tsne

setMethod("run.tsne", "spCounts", function(
    x,
    plot.callback,
    k = 3,
    max_iter = 5000,
    perplexity = 30,
    initial_dims = 50,
    ...
){
    counts.log <- getData(counts.log)
    maxs <- order(apply(counts.log, 1, max), decreasing=T)
    my.dist <- as.dist(1-cor(2^counts.log[maxs[1:2000], ], method="p"))
    my.tsne <- tsne(
        my.dist,
        k = k,
        epoch_callback = plot.callback,
        initial_dims = initial_dims,
        max_iter = max_iter,
        perplexity = perplexity,
        ...
    )
    rownames(my.tsne) <- rownames(my.dist)
})


