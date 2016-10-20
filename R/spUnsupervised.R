
#'@include All-classes.R
NULL

#' spUnsupervised
#'
#' Subtitle
#'
#' Description.
#'
#' @name spUnsupervised
#' @rdname spUnsupervised
#' @aliases spUnsupervised
#' @param spCounts spCounts object.
#' @param theta Passed to Rtsne.
#' @param k The dimensions of the resulting tsne. Passed to tsne function.
#' @param max_iter The max number of tSNE iterations. Passed to tsne function.
#' @param perplexity The perplexity argument to tsne. Passed to tsne function.
#' @param initial_dims The initial dimensions argument. Passed to tsne function.
#' @param Gmax A numeric vector of 1:Gmax passed as the "G" argument to Mclust.
#' @param seed Sets the seed before running tSNE.
#' @param type Decides if genes included are picked by their maximum expression
#'  or maximum variance. Can be either "max", "var", or "manual". If "manual" the
#'  genes argument must also be specified.
#' @param max The max number of genes to include based either on maximum expression
#'    or maximum variance which is decided by the "type" paramater.
#' @param genes If type = manual, genes to be included are specified here as a character
#'    vector. These must match the rownames in the counts variable.
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
#' @importFrom Rtsne Rtsne
#' @importFrom mclust Mclust mclustBIC
#' @importFrom plyr ddply summarize
#' @importFrom stats as.dist

setMethod("spUnsupervised", "spCounts", function(
    spCounts,
    theta = 0,
    k = 2,
    max_iter = 20000,
    perplexity = 10,
    initial_dims = 50,
    Gmax = 50,
    seed = 11,
    type = "max",
    max = 2000,
    genes=NULL,
    ...
){
    #get relevant data
    counts.log <- getData(spCounts, "counts.log")
    sampleType <- getData(spCounts, "sampleType")
    counts <- getData(spCounts, "counts")
    
    #filter genes to be included in analysis
    if(type == "var") {
        select <- .ntopVar(counts.log[ ,sampleType=="Singlet"], max)
    }
    
    if(type == "max") {
        select <- .ntopMax(counts.log[ ,sampleType=="Singlet"], max)
    }
    
    if(type == "manual" & is.character(genes) == TRUE) {
        select <- which(rownames(counts) %in% genes)
    }
    
    #calculate distances
    my.dist <- .distFunc(counts.log, select, sampleType)
    
    #run tSNE
    my.tsne <- .runTsne(my.dist, k, theta, initial_dims, max_iter, perplexity, seed)
    
    #run Mclust
    mod1 <- .runMclust(seed, my.tsne, Gmax)
    
    #create object
    new("spUnsupervised",
        counts=counts,
        counts.log=counts.log,
        sampleType=sampleType,
        tsne=my.tsne,
        tsneMeans=.tsneGroupMeans(my.tsne, mod1$classification),
        groupMeans=.averageGroupExpression(
            mod1$classification,
            counts[, sampleType == "Singlet"]
        ),
        classification=mod1$classification,
        selectInd=select
    )
})

#calculate the genes with the max variance and return the index
#for the top genes according to the value of "n".
.ntopVar <- function(data, n) {
    rv = apply(data, 1, var)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
}

#calculate the genes with the max expression and return the index
#for the top genes according to the value of "n".
.ntopMax <- function(data, n) {
    rv = apply(data, 1, max)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
}

#calculate the average expression of each gene within each classification group.
.averageGroupExpression <- function(classes, sng) {
    c <- unique(classes)
    means <- lapply(c, function(x) {
        rowMeans(sng[,c == x])
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
}

#calculate the mean x and y value for each classification group based
#on the tSNE results. Subsequently used for plotting.
.tsneGroupMeans <- function(x, class) {
    d <- data.frame(x[ ,1], x[ ,2], class)
    colnames(d) <- c("x", "y", "classification")
    means <- ddply(d, "classification", summarize, meanX=mean(x), meanY=mean(y))
    return(means)
}

#calculates the distance of 1-correlation
.distFunc <- function(x, select, sampleType) {
    as.dist(
        1-cor(
            2^x[select, sampleType == "Singlet"],
            method="p"
        )
    )
}

#runs the tSNE function
.runTsne <- function(my.dist, k, theta, initial_dims, max_iter, perplexity, seed) {
    set.seed(seed)
    
    my.tsne <- Rtsne(
        my.dist,
        k = k,
        #epoch_callback = plot.callback,
        initial_dims = initial_dims,
        max_iter = max_iter,
        perplexity = perplexity,
        theta = theta
    )$Y
    
    rownames(my.tsne) <- rownames(my.dist)
    return(my.tsne)
}

#runs the Mclust function
.runMclust <- function(seed, my.tsne, Gmax) {
    set.seed(seed)
    mod1 <- Mclust(my.tsne, G=1:Gmax)
    
    #rename classification classes
    x <- unique(mod1$classification)
    n <- ceiling(length(x)/26)
    names <- unlist(lapply(1:n, function(u) paste(LETTERS, u, sep="")))[1:length(x)]
    mod1$classification <- names[match(mod1$classification, x)]
    
    return(mod1)
}




