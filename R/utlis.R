
#'@include All-classes.R
NULL

#' selectTopVar
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' variance in the spCounts object. The expression matrix in
#' the counts.cpm slot is used for the calculation.
#'
#' @name selectTopVar
#' @rdname selectTopVar
#' @param cpm matrix; Counts per million.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @examples
#' s <- selectTopVar(getData(CIMseqSinglets_test, "counts.cpm"), 10)
NULL

#' @rdname selectTopVar
#' @importFrom matrixStats rowVars
#' @export

selectTopVar <- function(cpm, n) {
  rv <- matrixStats::rowVars(cpm)
  select <- order(rv, decreasing = TRUE)[1:n]
  return(select)
}

#' selectTopMax
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' expression in the spCounts object. The expression matrix in
#' the counts.cpm slot is used for the calculation.
#'
#' @name selectTopMax
#' @rdname selectTopMax
#' @param cpm matrix; Counts per million.
#' @param n Number of genes to select. If n > dim(data)[1] all data is returned.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @examples
#' s <- selectTopMax(getData(CIMseqSinglets_test, "counts.cpm"), 10)
NULL

#' @rdname selectTopMax
#' @importFrom matrixStats rowMaxs
#' @export

selectTopMax <- function(cpm, n) {
  n <- min(n, dim(cpm)[1])
  rv <- matrixStats::rowMaxs(cpm)
  select <- order(rv, decreasing = TRUE)[1:n]
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
#' @param cpm matrix; Counts per million.
#' @param select A numeric vector indicating the indexes of genes to include.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords pearsonsDist
#' @examples
#'
#' d <- pearsonsDist(getData(CIMseqSinglets_test, "counts.cpm"), 1:2000)
#'
NULL

#' @rdname pearsonsDist
#' @importFrom stats cor as.dist
#' @export

pearsonsDist <- function(cpm, select) {
  as.dist(
    1 - cor(
      cpm[select, ],
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
#' my.dist <- pearsonsDist(getData(CIMseqSinglets_test, "counts.cpm"), 1:2000)
#' tsne <- runTsne(my.dist, max_iter = 10)
#'
NULL

#' @rdname runTsne
#' @importFrom Rtsne Rtsne
#' @export

runTsne <- function(
  my.dist, dims = 2, theta = 0, initial_dims = 50, max_iter = 2000,
  perplexity = 10, seed = 11, is_distance = TRUE, ...
){
  set.seed(seed)
  
  my.tsne <- Rtsne(
    my.dist, dims = dims, initial_dims = initial_dims, max_iter = max_iter,
    perplexity = perplexity, theta = theta, is_distance = is_distance
  )$Y
  
  rownames(my.tsne) <- attr(my.dist, "Labels")
  return(my.tsne)
}

#' means.dim.red
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group in 
#' the dimensionality reduced data.
#'
#'
#' This method is typically only used in conjunction with plotting. It
#' calculates the 2 dimensional location of the mean of each classified group
#' in the supplied unsupervised dimensionality reduction (t-SNE) data
#' representation.
#'
#' @name means.dim.red
#' @rdname means.dim.red
#' @param data Singlet 2D tsne.
#' @param classes A character vector indicating the class of each singlet.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @examples
#'
#' means <- means.dim.red(
#'   getData(CIMseqSinglets_test, "dim.red"),
#'   getData(CIMseqSinglets_test, "classification")
#' )
#'
NULL

#' @rdname means.dim.red
#' @importFrom dplyr "%>%" group_by summarise
#' @importFrom rlang .data
#' @export

means.dim.red <- function(data, classes) {
  d <- data.frame(data[ ,1], data[ ,2], classes)
  colnames(d) <- c("x", "y", "classification")
  
  d %>%
  group_by(.data$classification) %>%
  summarise(
    x = mean(.data$x),
    y = mean(.data$y)
  ) %>%
  as.data.frame()
}

#' erccPerClass
#'
#' Calculates median and mean ercc reads per class.
#'
#' @name erccPerClass
#' @rdname erccPerClass
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @return A matrix containing the median and mean ERCC fractions for each
#'  classification group.
#' @author Jason T. Serviss
#' @examples
#'
#' out <- erccPerClass(CIMseqSinglets_test, CIMseqMultiplets_test)
#'
NULL

#' @rdname erccPerClass
#' @importFrom dplyr group_by summarise right_join
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export

erccPerClass <- function(
  singlets, multiplets
){
  estimateCells(singlets, multiplets) %>%
    right_join(tibble(
      class = getData(singlets, "classification"),
      sample = colnames(getData(singlets, "counts"))
    ), by = "sample") %>%
    group_by(class) %>%
    summarise(
      medianFracErcc = median(.data$frac.ercc, na.rm = TRUE),
      meanFracErcc = mean(.data$frac.ercc, na.rm = TRUE)
    )
}

#' cellNumberPerClass
#'
#' Calculates median and mean cell estimated cell number per class.
#'
#' @name cellNumberPerClass
#' @rdname cellNumberPerClass
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @return A matrix containing the median and mean cell numbers (estimated via
#' the ERCC fractions) for each classification group.
#' @author Jason T. Serviss
#' @examples
#'
#' out <- erccPerClass(CIMseqSinglets_test, CIMseqMultiplets_test)
#'
NULL

#' @rdname cellNumberPerClass
#' @importFrom dplyr group_by summarise right_join
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export

cellNumberPerClass <- function(
  singlets, multiplets
){
  estimatedCellNumber <- NULL
  estimateCells(singlets, multiplets) %>%
    right_join(tibble(
      class = getData(singlets, "classification"),
      sample = colnames(getData(singlets, "counts"))
    ), by = "sample") %>%
    group_by(class) %>%
    summarise(
      medianCellNumber = median(estimatedCellNumber),
      meanCellNumber = mean(estimatedCellNumber)
    )
}

#' countsPerClass
#'
#' Calculates median counts per class.
#'
#' @name countsPerClass
#' @rdname countsPerClass
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @return A matrix containing the median value for each gene for each
#'  classification group.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- countsPerClass(CIMseqSinglets_test)
#'
NULL

#' @rdname countsPerClass
#' @importFrom dplyr group_by summarise right_join
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export

countsPerClass <- function(
  singlets
){
  sampleName <- NULL
  counts <- getData(singlets, "counts")
  tibble(
    class = getData(singlets, "classification"),
    sampleName = colnames(getData(singlets, "counts"))
  ) %>%
  group_by(class) %>%
  summarize(
    medianCounts = median(colSums(counts[, colnames(counts) %in% sampleName]))
  )
}

#' expectedInteractionFreq
#'
#' Calculates the expected interaction frequency between cell types.
#'
#' @name expectedInteractionFreq
#' @rdname expectedInteractionFreq
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @return A table with the expected frequencies.
#' @author Jason T. Serviss
#' @keywords expectedInteractionFreq
#' @examples
#'
#' exfreq <- expectedInteractionFreq(CIMseqSinglets_test)
#'
NULL

#' @rdname expectedInteractionFreq
#' @importFrom dplyr "%>%"
#' @export

expectedInteractionFreq <- function(singlets) {
  getData(singlets, "classification") %>%
    table() %>% 
    divide_by(sum(.))
}

#' estimateTotalConnections
#'
#' Estimates the total number of connections in the multiplets.
#'
#' @name estimateTotalConnections
#' @rdname estimateTotalConnections
#' @aliases estimateTotalConnections
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @return The number of estimated connections.
#' @author Jason T. Serviss
#' @keywords estimateTotalConnections
#' @examples
#'
#' estimateTotalConnections(CIMseqSinglets_test, CIMseqMultiplets_test)
#'
NULL

#' @rdname estimateTotalConnections
#' @importFrom dplyr "%>%"
#' @importFrom utils combn
#' @export

estimateTotalConnections <- function(singlets, multiplets) {
  sampleType <- cellNumberMedian <- cellNumber <- connections <- NULL
  estimatedCellNumber <- NULL
  
  estimateCells(singlets, multiplets) %>%
    filter(sampleType == "Multiplet") %>%
    mutate(cellNumber = round(estimatedCellNumber)) %>%
    filter(cellNumber > 1) %>%
    mutate(connections = map_dbl(cellNumber, function(n) {
      ncol(combn(1:n, 2))
    })) %>%
    pull(connections) %>%
    sum()
}