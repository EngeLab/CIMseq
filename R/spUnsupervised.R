
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
#' @param spCounts spCounts object with singlets only.
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
#' @param counts Passed from spCounts object.
#' @param counts.log Passed from spCounts object.
#' @param sampleType Passed from spCounts object.
#' @param tsne tSNE results.
#' @param tsneMeans The mean x and y positions of each cell type in the tSNE results.
#' @param groupMeans The mean gene expression values for each cell type (classificaiton).
#' @param classification Post-tSNE cell type classification. Typically determined by the mclust package.
#' @param selectInd The indexes of the genes picked for use in spUnsupervised. Used in spSwarm.
#' @param object spUnsupervised object.
#' @param n Data to extract from spUnsupervised object.
#' @param value Data to replace in spUnsupervised object.
#' @param .Object Internal object.
#' @param ... Additional arguments to pass on
#' @return spUnsupervised object.
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
    max_iter = 2000,
    perplexity = 10,
    initial_dims = 50,
    Gmax = 50,
    seed = 11,
    type = "max",
    max = 2000,
    genes=NULL,
    ...
){
    
    #filter genes to be included in analysis
    if(type == "var") {
        select <- spTopVar(spCounts, max)
    }
    
    if(type == "max") {
        select <- spTopMax(spCounts, max)
    }
    
    #if(type == "manual" & is.character(genes) == TRUE) {
    #    select <- which(rownames(counts) %in% genes)
    #}
    
    #calculate distances
    my.dist <- pearsonsDist(spCounts, select)
    
    #run tSNE
    tsne <- runTsne(my.dist, k, theta, initial_dims, max_iter, perplexity, seed)
    
    #run Mclust
    class <- runMclust(tsne, Gmax, seed)
    
    #check for classification problems
    .classificationChecks(class)
    
    #create object
    new("spUnsupervised",
        tsne=tsne,
        tsneMeans=tsneGroupMeans(tsne, class),
        groupMeans=averageGroupExpression(spCounts,class),
        classification=class,
        selectInd=select
    )
})

#checks for results from classification that will throw a downstream error.
.classificationChecks <- function(class) {
    if(length(unique(class)) == 1) {
        stop("Only one group could be classified.
        Please adjust the tSNE-related arguments and try again or manually supply group classifications.")
    }
    
    if(any(table(class) == 1)) {
        stop("One/more cells have been classified as an individual cell type. 
        You probably need to adjust the Gmax argument and run again.")
    }
}



#' spTopVar
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' variance in the spCounts object. The expression matrix in
#' the counts.cpm slot is used for the calculation.
#'
#' @name spTopVar
#' @rdname spTopVar
#' @aliases spTopVar
#' @param spCounts An spCounts object.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords spTopVar
#' @examples
#'
#' cObj <- spCounts(testCounts, testErcc)
#' selected <- spTopVar(cObj, 10)
#'
NULL

#' @rdname spTopVar
#' @export


spTopVar <- function(spCounts, n) {
    rv = apply(getData(spCounts, "counts.cpm"), 1, var)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
}

#' spTopMax
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' expression in the spCounts object. The expression matrix in
#' the counts.cpm slot is used for the calculation.
#'
#' @name spTopMax
#' @rdname spTopMax
#' @aliases spTopMax
#' @param spCounts An spCounts object.
#' @param n Number of genes to select. If n > dim(data)[1] all data is returned.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords spTopMax
#' @examples
#'
#' cObj <- spCounts(testCounts, testErcc)
#' selected <- spTopMax(cObj, 10)
#'
NULL

#' @rdname spTopMax
#' @export

spTopMax <- function(spCounts, n) {
    data <- getData(spCounts, "counts.cpm")
    n <- min(n, dim(data)[1])
    rv = apply(data, 1, max)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
}

#' pearsonsDist
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group.
#'
#'
#' This method is typically only used in conjunction with plotting. It
#' calculates the 2 dimensional location of the mean of each classified group
#' in the supplied unsupervised dimensionality reduction (t-SNE) data
#' representation.
#'
#' @name pearsonsDist
#' @rdname pearsonsDist
#' @aliases pearsonsDist
#' @param data An spCounts object.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each classification group.
#' @author Jason T. Serviss
#' @keywords pearsonsDist
#' @examples
#'
#' singlets <- testCounts[ , grepl("^s.*", colnames(testCounts))]
#' select <- spTopMax(singlets, 10)
#' my.dist <- pearsonsDist(singlets, select)
#'
NULL

#' @rdname pearsonsDist
#' @export

pearsonsDist <- function(data, select) {
    as.dist(
        1-cor(
            2^getData(data, "counts.log")[select, ],
            method="p"
        )
    )
}

#' runTsne
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group.
#'
#'
#' This method is typically only used in conjunction with plotting. It
#' calculates the 2 dimensional location of the mean of each classified group
#' in the supplied unsupervised dimensionality reduction (t-SNE) data
#' representation.
#'
#' @name runTsne
#' @rdname runTsne
#' @aliases runTsne
#' @param data Singlet expression matrix.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each classification group.
#' @author Jason T. Serviss
#' @keywords runTsne
#' @examples
#'
#'\dontrun{
#' singlets <- testCounts[ , grepl("^s.*", colnames(testCounts))]
#' select <- spTopMax(singlets, 10)
#' my.dist <- pearsonsDist(singlets, select)
#' tsne <- runTsne(my.dist, max_iter=10)
#'}
#'
NULL

#' @rdname runTsne
#' @importFrom Rtsne Rtsne
#' @export

runTsne <- function(
    my.dist,
    k = 2,
    theta = 0,
    initial_dims = 50,
    max_iter = 2000,
    perplexity = 10,
    seed = 11
){
    set.seed(seed)
    
    my.tsne <- Rtsne(
        my.dist,
        k = k,
        initial_dims = initial_dims,
        max_iter = max_iter,
        perplexity = perplexity,
        theta = theta
    )$Y
    
    rownames(my.tsne) <- rownames(my.dist)
    return(my.tsne)
}

#' runMclust
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group.
#'
#'
#' This method is typically only used in conjunction with plotting. It
#' calculates the 2 dimensional location of the mean of each classified group
#' in the supplied unsupervised dimensionality reduction (t-SNE) data
#' representation.
#'
#' @name runMclust
#' @rdname runMclust
#' @aliases runMclust
#' @param data Singlet expression matrix.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each classification group.
#' @author Jason T. Serviss
#' @keywords runMclust
#' @examples
#'
#'\dontrun{
#' singlets <- testCounts[ , grepl("^s.*", colnames(testCounts))]
#' select <- spTopMax(singlets, 10)
#' my.dist <- pearsonsDist(singlets, select)
#' tsne <- runTsne(my.dist, max_iter=10)
#' class <- runMclust(tsne, 50, 11)
#'}
#'
NULL

#' @rdname runTsne
#' @importFrom mclust Mclust mclustBIC
#' @export

runMclust <- function(
    my.tsne,
    Gmax = 50,
    seed = 11
){
    set.seed(seed)
    mod1 <- Mclust(my.tsne, G=1:Gmax)
    
    #rename classification classes
    x <- unique(mod1$classification)
    n <- ceiling(length(x)/26)
    names <- unlist(lapply(1:n, function(u) paste(LETTERS, u, sep="")))[1:length(x)]
    mod1$classification <- names[match(mod1$classification, x)]
    
    return(mod1$classification)
}

#' averageGroupExpression
#'
#' This output from this function is utilized to represent each group classified
#' group during swarm optimization.
#'
#' Calculate the average expression of each gene within each classification group.
#' Note that typically only singlets are input to this method.
#'
#' @name averageGroupExpression
#' @rdname averageGroupExpression
#' @aliases averageGroupExpression
#' @param data Singlet expression matrix.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each classification group.
#' @author Jason T. Serviss
#' @keywords averageGroupExpression
#' @examples
#'
#'\dontrun{
#' singletIdx <- grepl("^s.*", colnames(testCounts))
#' sngMatrix <- testCounts[ , singletIdx]
#' sngErcc <- testErcc[ , singletIdx]
#' cObj <- spCounts(sngMatrix, sngErcc)
#' select <- spTopMax(cObj, 10)
#' my.dist <- pearsonsDist(cObj, select)
#' tsne <- runTsne(my.dist, max_iter=2000)
#' class <- runMclust(tsne, 50, 11)
#' averageExp <- averageGroupExpression(cObj, class)
#'}
#'
NULL

#' @rdname averageGroupExpression
#' @export

averageGroupExpression <- function(data, classes) {
    c <- unique(classes)
    exp <- getData(data, "counts.cpm")
    means <- lapply(c, function(x) {
        rowMeans(exp[ ,classes == x])
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
}

#' tsneGroupMeans
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group.
#'
#'
#' This method is typically only used in conjunction with plotting. It
#' calculates the 2 dimensional location of the mean of each classified group
#' in the supplied unsupervised dimensionality reduction (t-SNE) data
#' representation.
#'
#' @name tsneGroupMeans
#' @rdname tsneGroupMeans
#' @aliases tsneGroupMeans
#' @param data Singlet expression matrix.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each classification group.
#' @author Jason T. Serviss
#' @keywords tsneGroupMeans
#' @examples
#'
#' singlets <- testCounts[ , grepl("^s.*", colnames(testCounts))]
#' tsne <- runTsne(singlets, max_iter=10)
#' classes <- runMclust(tsne, 50, 11)
#' meanGroupPos <- tsneGroupMeans(singlets, classes)
#'
NULL

#' @rdname tsneGroupMeans
#' @importFrom plyr ddply
#' @export

tsneGroupMeans <- function(data, classes) {
    d <- data.frame(data[ ,1], data[ ,2], classes)
    colnames(d) <- c("x", "y", "classification")
    means <- ddply(
        d,
        "classification",
        summarize,
        meanX=mean(substitute(x)),
        meanY=mean(substitute(y))
    )
    return(means)
}






