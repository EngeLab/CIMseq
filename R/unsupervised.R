
#'@include All-classes.R
NULL

#' Run Tsne
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name spUnsupervised
#' @rdname spUnsupervised
#' @aliases spUnsupervised
#' @param counts Counts object.
#' @param ... Additional arguments to pass on
#' @return Ercc fraction plot.
#' @author Jason T. Serviss
#' @keywords spUnsupervised
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname spUnsupervised
#' @export

setGeneric("spUnsupervised", function(spCounts, ...
){ standardGeneric("spUnsupervised") })


#' @rdname spUnsupervised
#' @export
#' @importFrom tsne tsne
#' @importFrom mclust Mclust mclustBIC

setMethod("spUnsupervised", "spCounts", function(
    spCounts,
    max = 2000,
    plot.callback = NULL,
    k = 2,
    max_iter = 10000,
    perplexity = 10,
    initial_dims = 50,
    Gmax = 50,
    seed = 11,
    ...
){
    counts.log <- getData(spCounts, "counts.log")
    sampleType <- getData(spCounts, "sampleType")
    
    maxs <- order(apply(counts.log, 1, max), decreasing=T) ##use feature selection here instead of highest expressed?
    my.dist <- as.dist(
        1-cor(
            2^counts.log[maxs[1:max],
            sampleType != "Doublet"],
            method="p")
    )
    
    set.seed(seed)
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
    
    set.seed(seed)
    mod1 <- Mclust(my.tsne, G=1:Gmax)
    
    new("spUnsupervised",
        counts.log=counts.log,
        dist=as.matrix(my.dist),
        tsne=my.tsne,
        mclust=unclass(mod1)
    )
})

.ntopF <- function(data = y, n = ntop) {
    if(is.numeric(n)){
        rv = genefilter::rowVars(data)
        select = order(rv, decreasing=TRUE)[seq_len(min(n, length(rv)))]
        return(select)
    }
}
