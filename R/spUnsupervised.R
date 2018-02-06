
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
#'  or maximum variance. Can be either "max", "var", "maxMean", or "manual". If
#'  "manual" the genes argument must also be specified.
#' @param max The max number of genes to include based either on maximum
#'    expression or maximum variance which is decided by the "type" paramater.
#' @param genes If type = manual, genes to be included are specified here as a
#'    character vector. These must match the rownames in the counts variable.
#' @param tsne tSNE results.
#' @param tsneMeans The mean x and y positions of each cell type in the tSNE
#'    results.
#' @param groupMeans The mean gene expression values for each cell type
#'    (classificaiton).
#' @param classification Post-tSNE cell type classification. Typically
#'    determined by the mclust package.
#' @param selectInd The indexes of the genes picked for use in spUnsupervised.
#' @param weighted boolean indicating if the group means should be weighted with
#'    uncertainty.
#' @param uncertainty Contains the uncertainty in converting a conditional
#'    probablility from EM to a classification in model-based clustering.
#'    Reported from mclust.
#' @param object spUnsupervised object.
#' @param n Data to extract from spUnsupervised object.
#' @param value Data to replace in spUnsupervised object.
#' @param .Object Internal object.
#' @param ... Additional arguments to pass on.
#' @return spUnsupervised object.
#' @author Jason T. Serviss
#' @keywords spUnsupervised
#' @examples
#'
#' #use test data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#'
#' #run function
#' uObj <- spUnsupervised(cObjSng, max_iter = 100, Gmax = 6)
#'
#' #run function
#'
NULL

#' @rdname spUnsupervised
#' @export

setGeneric("spUnsupervised", function(
    spCounts,
    ...
){
    standardGeneric("spUnsupervised")
})


#' @rdname spUnsupervised
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom mclust Mclust mclustBIC
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
    genes = NULL,
    weighted = TRUE,
    ...
){
    #filter genes to be included in analysis
    select <- .featureSelection(spCounts, type, max, genes)
    
    #calculate distances
    my.dist <- pearsonsDist(spCounts, select)
    
    #run tSNE
    tsne <- runTsne(
        my.dist,
        k,
        theta,
        initial_dims,
        max_iter,
        perplexity,
        seed
    )
    
    #run Mclust
    tmp <- runMclust(tsne, Gmax, seed)
    class <- tmp[[1]]
    uncertainty <- tmp[[2]]
    
    #check for classification problems
    .classificationChecks(class)
    
    #create object
    new("spUnsupervised",
        tsne = tsne,
        tsneMeans = tsneGroupMeans(tsne, class),
        groupMeans = averageGroupExpression(
            spCounts,
            class,
            weighted,
            uncertainty
        ),
        classification = class,
        uncertainty = uncertainty,
        selectInd = select
    )
})

.featureSelection <- function(spCounts, type, max, genes) {
  if(type == "var") {
    select <- spTopVar(spCounts, max)
  }
  
  if(type == "max") {
    select <- spTopMax(spCounts, max)
  }
    
  if(type == "maxMean") {
    select <- spTopMaxMean(spCounts, max)
  }
    
  if(type == "manual" & is.character(genes) == TRUE) {
    select <- which(rownames(getData(spCounts, "counts.log")) %in% genes)
  }
    
  if(type == "all") {
    bool1 <- duplicated(getData(spCounts, "counts.log"))
    bool2 <- duplicated(getData(spCounts, "counts.log"), fromLast = TRUE)
    select <- which((bool1 | bool2) == FALSE)
  }
  return(select)
}
#checks for results from classification that will throw a downstream error.
.classificationChecks <- function(class) {
    if(length(unique(class)) == 1) {
        stop("Only one group could be classified.
        Please adjust the tSNE-related arguments and
        try again or manually supply group classifications.")
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
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' selected <- spTopVar(cObjSng, 10)
#'
NULL

#' @rdname spTopVar
#' @importFrom stats var
#' @export

spTopVar <- function(spCounts, n) {
    rv = apply(getData(spCounts, "counts.cpm"), 1, var)
    select = order(rv, decreasing = TRUE)[1:n]
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
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' selected <- spTopMax(cObjSng, 10)
#'
NULL

#' @rdname spTopMax
#' @export

spTopMax <- function(spCounts, n) {
    data <- getData(spCounts, "counts.cpm")
    n <- min(n, dim(data)[1])
    rv <- apply(data, 1, max)
    select <- order(rv, decreasing = TRUE)[1:n]
    return(select)
}

#' spTopMaxMean
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' max(gene) - mean(gene) values in the spCounts object. The expression matrix
#' in the counts.cpm slot is used for the calculation.
#'
#' @name spTopMaxMean
#' @rdname spTopMaxMean
#' @aliases spTopMaxMean
#' @param spCounts An spCounts object.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords spTopMaxMean
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' selected <- spTopMaxMean(cObjSng, 10)
#'
NULL

#' @rdname spTopVar
#' @importFrom stats var
#' @export

spTopMaxMean <- function(spCounts, n) {
  rv = apply(getData(spCounts, "counts.cpm"), 1, function(x) max(x) - mean(x))
  select = order(rv, decreasing = TRUE)[1:n]
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
#' @param spCounts An spCounts object.
#' @param select A numeric vector indicating the indexes of genes to include.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords pearsonsDist
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' my.dist <- pearsonsDist(cObjSng, 1:nrow(testCounts))
#'
NULL

#' @rdname pearsonsDist
#' @importFrom stats cor
#' @export

pearsonsDist <- function(spCounts, select) {
    as.dist(
        1-cor(
            getData(spCounts, "counts.log")[select, ],
            method = "p"
        )
    )
}

#' runTsne
#'
#'
#' Calculates the x and y coordinates of the mean of each classif	ied group.
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
#' @param my.dist A distance object typically produced with pearsonsDist.
#' @param dims Argument to Rtsne. Numeric indicating the output dimensions.
#' @param theta Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param initial_dims Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param max_iter Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param perplexity Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param seed The desired seed to set before running.
#' @param is_distance Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param ... Additional arguments to pass on
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords runTsne
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[ ,s], testErcc[, s])
#' select <- spTopMax(cObjSng, 10)
#' my.dist <- pearsonsDist(cObjSng, select)
#' tsne <- runTsne(my.dist, max_iter = 10)
#'
#'
NULL

#' @rdname runTsne
#' @importFrom Rtsne Rtsne
#' @export

runTsne <- function(
    my.dist,
    dims = 2,
    theta = 0,
    initial_dims = 50,
    max_iter = 2000,
    perplexity = 10,
    seed = 11,
    is_distance = TRUE,
    ...
){
    set.seed(seed)
    
    my.tsne <- Rtsne(
        my.dist,
        dims = dims,
        initial_dims = initial_dims,
        max_iter = max_iter,
        perplexity = perplexity,
        theta = theta,
        is_distance = is_distance
    )$Y
    
    rownames(my.tsne) <- attr(my.dist, "Labels")
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
#' @param my.tsne A tsne result.
#' @param Gmax A integer indicating the maximum number of clusters to evaluate.
#' @param seed The desired seed to set before running.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords runMclust
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' select <- spTopMax(cObjSng, 10)
#' my.dist <- pearsonsDist(cObjSng, select)
#' tsne <- runTsne(my.dist, max_iter = 10)
#' class <- runMclust(tsne, 6, 11)
#'
NULL

#' @rdname runMclust
#' @importFrom mclust Mclust mclustBIC
#' @export

runMclust <- function(
    my.tsne,
    Gmax = 50,
    seed = 11
){
    set.seed(seed)
    mod1 <- Mclust(my.tsne, G = 1:Gmax)
    
    #rename classification classes
    x <- unique(mod1$classification)
    n <- ceiling(length(x)/26)
    names <- unlist(lapply(1:n, function(u)
        paste(LETTERS, u, sep = ""))
    )[1:length(x)]
    mod1$classification <- names[match(mod1$classification, x)]
    
    classification <- mod1$classification
    uncertainty <- mod1$uncertainty
    return(list(classification, uncertainty))
}

#' averageGroupExpression
#'
#' This output from this function is utilized to represent each group classified
#' group during swarm optimization.
#'
#' Calculate the average expression of each gene within each classification
#' group. Note that typically only singlets are input to this method.
#'
#' @name averageGroupExpression
#' @rdname averageGroupExpression
#' @aliases averageGroupExpression
#' @param data Singlet expression matrix.
#' @param classes A character vector indicating the class of each singlet.
#' @param weighted Logical indicating if the group means shoule be weighted with
#'    the uncertainty.
#' @param uncertainty A numeric vector indicating the uncertainty if weighted is
#'    TRUE.
#' @return A matrix containing the mean value for each gene for each
#' classification group.
#' @author Jason T. Serviss
#' @keywords averageGroupExpression
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' select <- spTopMax(cObjSng, 10)
#' my.dist <- pearsonsDist(cObjSng, select)
#' tsne <- runTsne(my.dist, max_iter = 2000)
#' mclustOUT <- runMclust(tsne, 7, 11)
#' class <- mclustOUT[[1]]
#' uncertainty <- mclustOUT[[2]]
#'
#' #Weighted mean
#' weighted <- TRUE
#' averageExp <- averageGroupExpression(cObjSng, class, weighted, uncertainty)
#'
#' #Unweighted mean
#' weighted <- FALSE
#' averageExp <- averageGroupExpression(cObjSng, class, weighted)
#'
NULL

#' @rdname averageGroupExpression
#' @export

averageGroupExpression <- function(
    data,
    classes,
    weighted,
    uncertainty = NULL
){
    c <- unique(classes)
    exp <- getData(data, "counts.cpm")
    
    if(weighted) {
        u <- 1 - uncertainty
        w <- t(t(exp) * u)
        means <- lapply(c, function(x) {
            rowSums(w[, classes == x]) / sum(u[classes == x])
        })
    } else {
        means <- lapply(c, function(x) {
            rowMeans(exp[, classes == x])
        })
    }
    
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
#' @param data Singlet 2D tsne.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords tsneGroupMeans
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' my.dist <- pearsonsDist(cObjSng)
#' tsne <- runTsne(my.dist, max_iter = 10)
#' classes <- runMclust(tsne, 7, 11)[[1]]
#' meanGroupPos <- tsneGroupMeans(tsne, classes)
#'
NULL

#' @rdname tsneGroupMeans
#' @import dplyr
#' @export

tsneGroupMeans <- function(data, classes) {
    d <- data.frame(data[ ,1], data[ ,2], classes)
    colnames(d) <- c("x", "y", "classification")
    
    means <- d %>%
        group_by(classification) %>%
        summarise(
            x = mean(.data$x),
            y = mean(.data$y)
        ) %>%
        as.data.frame()
        
    return(means)
}

#' erccPerClass
#'
#'
#' Calculates median ercc reads per class.
#'
#'
#'
#' @name erccPerClass
#' @rdname erccPerClass
#' @aliases erccPerClass
#' @param spCountsSng An spCounts object containing singlets.
#' @param spCountsMul An spCounts object containing multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords erccPerClass
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#' uObj <- testUns
#' output <- erccPerClass(cObjSng, cObjMul, uObj)
NULL

#' @rdname erccPerClass
#' @importFrom dplyr group_by summarise right_join
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export

erccPerClass <- function(
    spCountsSng,
    spCountsMul,
    spUnsupervised
){
    d <- estimateCells(spCountsSng, spCountsMul)
    class <- tibble(
        class = getData(spUnsupervised, "classification"),
        sampleName = colnames(getData(spCountsSng, "counts"))
    )
    
    d %>%
        right_join(class, by = "sampleName") %>%
        group_by(class) %>%
        summarise(medianFracErcc = median(.data$frac.ercc))
}





