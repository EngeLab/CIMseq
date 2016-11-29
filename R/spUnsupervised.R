
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
#' @param dims The dimensions of the resulting tsne. Passed to tsne function.
#' @param max_iter The max number of tSNE iterations. Passed to tsne function.
#' @param perplexity The perplexity argument to tsne. Passed to tsne function.
#' @param initial_dims The initial dimensions argument. Passed to tsne function.
#' @param Gmax A numeric vector of 1:Gmax passed as the "G" argument to Mclust.
#' @param seed Sets the seed before running tSNE.
#' @param geneSelectType Decides if genes included are picked by their maximum expression
#'  or maximum variance. Can be either "max", "var", or "manual". If "manual" the
#'  genes argument must also be specified.
#' @param max The max number of genes to include based either on maximum expression
#'    or maximum variance which is decided by the "type" paramater.
#' @param genes If type = manual, genes to be included are specified here as a character
#'    vector. These must match the rownames in the counts variable.
#' @param unsupervisedC tSNE results.
#' @param groupMeans The mean gene expression values for each cell type (classificaiton).
#' @param classification Post-tSNE cell type classification. Typically determined by the mclust package.
#' @param selectIdx The indexes of the genes picked for use in spUnsupervised. Used in spSwarm.
#' @param object spUnsupervised object.
#' @param n Data to extract from spUnsupervised object.
#' @param value Data to replace in spUnsupervised object.
#' @param x Default plot param, an spCounts object containing singlets.
#' @param y Default plot param, an spCounts object containing multuplets.
#' @param type Character; The type of plot desired. Currently \emph{markers} or \emph{ercc}.
#' @param markers Markers/genes to plot.
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

setMethod("spUnsupervised", "missing", function(
...
){
    new("spUnsupervised",
        unsupervisedC=matrix(),
        groupMeans=matrix(),
        classification=factor(),
        selectIdx=numeric()
    )
})

#' @rdname spUnsupervised
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom mclust Mclust mclustBIC
#' @importFrom plyr ddply summarize
#' @importFrom stats as.dist

setMethod("spUnsupervised", "spCounts", function(
    spCounts,
    theta = 0.0,
    dims = 2,
    max_iter = 20000,
    perplexity = 10,
    initial_dims = 50,
    Gmax = 50,
    seed = 11,
    geneSelectType = "max",
    max = 2000,
    genes=NULL,
    ...
){
    #get relevant data
    counts.log <- getData(spCounts, "counts.log")
    counts <- getData(spCounts, "counts")
    
    #filter genes to be included in analysis
    if(geneSelectType == "var") {
        select <- spNtopVar(counts.log, max)
    }
    
    if(geneSelectType == "max") {
        select <- spNtopMax(counts.log, max)
    }
    
    if(geneSelectType == "manual" & is.character(genes) == TRUE) {
        select <- which(rownames(counts) %in% genes)
    }
    
    #calculate distances
    my.dist <- spPearsonDist(counts.log, select)
    
    #run tSNE
    my.tsne <- spTsne(my.dist, dims, theta, initial_dims, max_iter, perplexity, seed)
    
    #run Mclust
    mod1 <- spMclust(my.tsne, seed, Gmax)
    
    #check for classification problems
    .classificationChecks(mod1)
    
    #create object
    new("spUnsupervised",
        unsupervisedC=my.tsne,
        groupMeans=spAverageGroupExpression(
            counts,
            mod1$classification
        ),
        classification=mod1$classification,
        selectIdx=select
    )
})

#checks for results from classification that will throw a downstream error.
.classificationChecks <- function(mod1) {
    if(length(unique(mod1$classification)) == 1) {
        stop("Only one group could be classified.
        Please adjust the tSNE-related arguments and try again or manually supply group classifications.")
    }
    
    if(any(table(mod1$classification) == 1)) {
        stop("One/more cells have been classified as an individual cell type. 
        You probably need to adjust the Gmax argument and run again.")
    }
}

#calculate the mean x and y value for each classification group based
#on the tSNE results. Subsequently used for plotting.
#tsneMeans=.tsneGroupMeans(my.tsne, mod1$classification),
##' @importFrom plyr ddply summarize

.tsneGroupMeans <- function(x, class) {
    d <- data.frame(x[ ,1], x[ ,2], class)
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


#' spNtopVar
#'
#' Subtitle
#'
#' Ranks genes by their variance and returns the index
#' for the number of genes with maximum variance according
#' to the value of \emph{n}.
#'
#' @name spNtopVar
#' @rdname spNtopVar
#' @aliases spNtopVar
#' @param counts A counts matrix or spCounts object.
#' @param n Integer; the number of gene indices to return.
#' @param ... Additional arguments to pass on
#' @return Numeric; Indices of selected genes.
#' @author Jason T. Serviss
#' @keywords spNtopVar
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname spNtopVar
#' @export

setGeneric("spNtopVar", function(counts, ...
){ standardGeneric("spNtopVar") })

#' @rdname spNtopVar
#' @export

setMethod("spNtopVar", "spCounts", function(
    counts,
    n,
    ...
){
    counts <- getData(counts, "counts.log")
    rv = apply(counts, 1, var)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
})

#' @rdname spNtopVar
#' @export

setMethod("spNtopVar", "matrix", function(
    counts,
    n,
    ...
){
    rv = apply(counts, 1, var)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
})

#' spNtopMax
#'
#' Subtitle
#'
#' Ranks genes by their expression and returns the index
#' for the number of genes with maximum expression according
#' to the value of \emph{n}.
#'
#' @name spNtopMax
#' @rdname spNtopMax
#' @aliases spNtopMax
#' @param counts A counts matrix or spCounts object.
#' @param n Integer; the number of gene indices to return.
#' @param ... Additional arguments to pass on
#' @return Numeric; Indices of selected genes.
#' @author Jason T. Serviss
#' @keywords spNtopMax
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname spNtopMax
#' @export

setGeneric("spNtopMax", function(counts, ...
){ standardGeneric("spNtopMax") })

#' @rdname spNtopMax
#' @export

setMethod("spNtopMax", "spCounts", function(
    counts,
    n,
    ...
){
    counts <- getData(counts, "counts.log")
    rv = apply(counts, 1, max)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
})

#' @rdname spNtopMax
#' @export

setMethod("spNtopMax", "matrix", function(
    counts,
    n,
    ...
){
    rv = apply(counts, 1, max)
    select = order(rv, decreasing=TRUE)[1:n]
    return(select)
})

#' spPearsonDist
#'
#' Subtitle
#'
#' Calculates 1-Pearson's correlation of samples using
#' all genes specified by \emph{selectIdx}.
#'
#' @name spPearsonDist
#' @rdname spPearsonDist
#' @aliases spPearsonDist
#' @param counts A counts matrix or spCounts object.
#' @param selectIdx Numeric; A vector indicating the indices of genes to use.
#' @param ... Additional arguments to pass on
#' @return Dist; containg 1-Pearson's correlation of samples.
#' @author Jason T. Serviss
#' @keywords spPearsonDist
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname spPearsonDist
#' @importFrom stats as.dist
#' @export

setGeneric("spPearsonDist", function(counts, ...
){ standardGeneric("spPearsonDist") })

#' @rdname spPearsonDist
#' @export

setMethod("spPearsonDist", "spCounts", function(
    counts,
    selectIdx,
    ...
){
    counts <- getData(counts, "counts.log")
    as.dist(
        1-cor(
            2^counts[selectIdx, ],
            method="p"
        )
    )
})

#' @rdname spPearsonDist
#' @importFrom stats as.dist
#' @export

setMethod("spPearsonDist", "matrix", function(
    counts,
    selectIdx,
    ...
){
    as.dist(
        1-cor(
            2^counts[selectIdx, ],
            method="p"
        )
    )
})

#' spTsne
#'
#' Runs Rtsne function from Rtsne package.
#'
#' This function can be used for unsupervised clustering of cells allowing
#' downstream cell classification based on the clusters formed.
#'
#' @name spTsne
#' @rdname spTsne
#' @aliases spTsne
#' @param distances Dist; potentially produced by \emph{spPearsonDist} function.
#' @param dims Integer; number of output dimensions. Passed to Rtsne.
#' @param theta Numeric; Speed/accuracy tradeoff. Value between 0 and 1, where
#' 0 is the most correct and 1 is the fastest. Passed to Rtsne.
#' @param initial_dims Integer; Number of dims to be retained in the initial
#' PCA step. Passed to Rtsne.
#' @param max_iter Integer; Number of iterations. Passed to Rtsne.
#' @param perplexity Numeric; Perplexity parameter. Passed to Rtsne.
#' @param seed Integer; Used to set seed.
#' @param is_distance Logical; Currently should not be altered. Passed to Rtsne.
#' @param ... Additional arguments to pass on
#' @return Matrix; new representations for the objects.
#' @author Jason T. Serviss
#' @keywords spTsne
#' @examples
#'
#' #use demo data
#' dist <- spPearsonDist(
#'      testCounts[ ,grepl("s.", colnames(testCounts))],
#'      selectIdx=1:nrow(testCounts)
#' )
#'
#' #run function
#' tsne <- spTsne(dist, max_iter=2, theta=1)
#'
NULL

#' @rdname spTsne
#' @export

setGeneric("spTsne", function(distances, ...
){ standardGeneric("spTsne") })

#' @rdname spTsne
#' @export
#' @importFrom Rtsne Rtsne

setMethod("spTsne", "dist", function(
    distances,
    dims = 2,
    theta = 0.0,
    initial_dims = 50,
    max_iter = 2000,
    perplexity = 10,
    seed = 11,
    is_distance=TRUE,
    ...
){
    set.seed(seed)
    
    my.tsne <- Rtsne(
        distances,
        dims = dims,
        initial_dims = initial_dims,
        max_iter = max_iter,
        perplexity = perplexity,
        theta = theta,
        is_distance=TRUE
    )$Y
    
    rownames(my.tsne) <- rownames(distances)
    return(my.tsne)
})

#' spAverageGroupExpression
#'
#' Calculates the average expression value of each gene in each class.
#'
#' This function calculates the mean expression of each gene in each class which
#' is necessary for downstream assignment of multuplets to specific singlet classes.
#'
#' @name spAverageGroupExpression
#' @rdname spAverageGroupExpression
#' @aliases spAverageGroupExpression
#' @param counts A counts matrix or spCounts object.
#' @param classes Character; A character vector indicating the class of each input sample.
#' @param ... Additional arguments to pass on
#' @return Matrix; Contains mean values for each gene in each class.
#' @author Jason T. Serviss
#' @keywords spAverageGroupExpression
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname spAverageGroupExpression
#' @export

setGeneric("spAverageGroupExpression", function(counts, ...
){ standardGeneric("spAverageGroupExpression") })

#' @rdname spAverageGroupExpression
#' @export

setMethod("spAverageGroupExpression", "spCounts", function(
    counts,
    classes,
    ...
){
    counts <- getData(counts, "counts")
    .spAveGrpExp.inputChecks(counts, classes)

    c <- unique(classes)
    means <- lapply(c, function(x) {
        rowMeans(counts[,classes == x])
    })
    
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
})

#' @rdname spAverageGroupExpression
#' @export

setMethod("spAverageGroupExpression", "matrix", function(
    counts,
    classes,
    ...
){
    .spAveGrpExp.inputChecks(counts, classes)
    
    c <- unique(classes)
    means <- lapply(c, function(x) {
        rowMeans(counts[,classes == x])
    })
    
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    return(means)
})

.spAveGrpExp.inputChecks <- function(counts, classes) {
    if(ncol(counts) != length(classes)) {
        stop("Several samples are missing a classification")
    }
}

#' spMclust
#'
#' Classification of unsupervised clustering data.
#'
#' This function calculates the mean expression of each gene in each class which
#' is necessary for downstream assignment of multuplets to specific singlet classes.
#'
#' @name spMclust
#' @rdname spMclust
#' @aliases spMclust
#' @param unsupervised Matrix; Contains matrix resulting from unsupervised clustering.
#' @param seed Integer; Used to set seed.
#' @param Gmax Integer; Indicates the maximum number of clusters for which BIC is
#' calculated. Passed to Mclust.
#' @param ... Additional arguments to pass on
#' @return Character; Contains classifications for each sample.
#' @author Jason T. Serviss
#' @keywords spMclust
#' @examples
#'
#' #use demo data
#'
#' #run function
#'
NULL

#' @rdname spMclust
#' @export

setGeneric("spMclust", function(unsupervised, ...
){ standardGeneric("spMclust") })

#' @rdname spMclust
#' @importFrom mclust Mclust mclustBIC
#' @export

setMethod("spMclust", "matrix", function(
    unsupervised,
    seed,
    Gmax,
    ...
){
    set.seed(seed)
    mod1 <- Mclust(unsupervised, G=1:Gmax)
    
    #rename classification classes
    x <- unique(mod1$classification)
    n <- ceiling(length(x)/26)
    names <- unlist(lapply(1:n, function(u) paste(LETTERS, u, sep="")))[1:length(x)]
    mod1$classification <- names[match(mod1$classification, x)]
    
    return(mod1)
})




