
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
#' @param pca matrix; Optional precomputed representation of the data in PCA
#' space.
#' @param pcs integer; Optional value of the largest principal component to
#' include from the PCA dimensionality reduction.
#' @param pcVarPercent numeric; Specifies the minimum amount of variance
#' contained in a principal component in order for it to be retained.
#' @param mask integer; Principal components to mask.
#' @param NN data.frame; Optional precalculated nearest neighbors. Should
#' contain columns "from" and "to" indicating neighbors, "dist" indicating a
#' distance between neighbors, and "mutual" a logical column indicating if the
#' neighbors have a mutual connection.
#' @param kNN integer; Number of nearest neighbors to identify.
#' @param distCut numeric; The quantile at which to cut the distance.
#' @param classCut integer; Classes with members < classCut will be "undefined".
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
#' uObj <- spUnsupervised(test_spCountsSng, max_iter = 100, distCut = 0.7, classCut = 1)
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
  spCounts, theta = 0, k = 2, max_iter = 2000, perplexity = 10,
  initial_dims = 50,
  pca = NULL,
  pcs = NULL,
  pcVarPercent = 1,
  mask = NULL,
  NN = NULL,
  kNN = 50,
  distCut,
  classCut,
  seed = 11, type = "max",
  max = 2000, genes = NULL, weighted = TRUE, ...
){
  ji <- NULL
  #filter genes to be included in analysis
  select <- .featureSelection(spCounts, type, max, genes)
  
  #calculate distances
  my.dist <- pearsonsDist(spCounts, select)
  
  #run tSNE
  tsne <- runTsne(
    my.dist, k, theta, initial_dims,
    max_iter, perplexity, seed
  )
  
  #run Mclust
  if(is.null(pca)) pca <- runPCA(spCounts, select)
  if(is.null(pcs)) pcs <- pcByVarPercent(pca, pcVarPercent)
  if(is.null(NN)) NN <- getKNN(pca, pcs, kNN, mask = mask)
  g <- constructGraph(NN, quantile(NN$dist, probs = distCut))
  g <- jaccardSimilarity(g, 0, seed)
  g <- classify_louvain(g, ji)
  g <- removeOutliers(g, "louvain", threshold = classCut)
  class <-  as.data.frame(g)$louvain
  
  #check for classification problems
  .classificationChecks(class)
  
  #create object
  new("spUnsupervised",
    tsne = tsne,
    tsneMeans = tsneGroupMeans(tsne, class),
    classification = class,
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
#' selected <- spTopVar(test_spCountsSng, 10)
#'
NULL

#' @rdname spTopVar
#' @importFrom matrixStats rowVars
#' @export

spTopVar <- function(spCounts, n) {
  rv <- matrixStats::rowVars(getData(spCounts, "counts.cpm"))
  select <- order(rv, decreasing = TRUE)[1:n]
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
#' selected <- spTopMax(test_spCountsSng, 10)
#'
NULL

#' @rdname spTopMax
#' @importFrom matrixStats rowMaxs
#' @export

spTopMax <- function(spCounts, n) {
  data <- getData(spCounts, "counts.cpm")
  n <- min(n, dim(data)[1])
  rv <- matrixStats::rowMaxs(data)
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
#' @param spCounts An spCounts object.
#' @param select A numeric vector indicating the indexes of genes to include.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords pearsonsDist
#' @examples
#'
#' my.dist <- pearsonsDist(test_spCountsSng, 1:2000)
#'
NULL

#' @rdname pearsonsDist
#' @importFrom stats cor
#' @export

pearsonsDist <- function(spCounts, select) {
  as.dist(
    1 - cor(
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
#' my.dist <- pearsonsDist(test_spCountsSng, 1:2000)
#' tsne <- runTsne(my.dist, max_iter = 10)
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
    my.dist, dims = dims, initial_dims = initial_dims, max_iter = max_iter,
    perplexity = perplexity, theta = theta, is_distance = is_distance
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
#' my.dist <- pearsonsDist(test_spCountsSng, 1:2000)
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
#' meanGroupPos <- tsneGroupMeans(
#'   getData(test_spUnsupervised, "tsne"),
#'   getData(test_spUnsupervised, "classification")
#' )
#'
NULL

#' @rdname tsneGroupMeans
#' @importFrom dplyr "%>%" group_by summarise
#' @export

tsneGroupMeans <- function(data, classes) {
  d <- data.frame(data[ ,1], data[ ,2], classes)
  colnames(d) <- c("x", "y", "classification")
  
  d %>%
  group_by(classification) %>%
  summarise(
    x = mean(.data$x),
    y = mean(.data$y)
  ) %>%
  as.data.frame()
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
#' @return A matrix containing the median value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords erccPerClass
#' @examples
#'
#' out <- erccPerClass(test_spCountsSng, test_spCountsMul, test_spUnsupervised)
#'
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
  class <- tibble(
    class = getData(spUnsupervised, "classification"),
    sampleName = colnames(getData(spCountsSng, "counts"))
  )
  
  estimateCells(spCountsSng, spCountsMul) %>%
  right_join(class, by = "sampleName") %>%
  group_by(class) %>%
  summarise(
    medianFracErcc = median(.data$frac.ercc, na.rm = TRUE),
    meanFracErcc = mean(.data$frac.ercc, na.rm = TRUE)
  )
}

#' countsPerClass
#'
#'
#' Calculates median counts per class.
#'
#'
#'
#' @name countsPerClass
#' @rdname countsPerClass
#' @aliases countsPerClass
#' @param spCountsSng An spCounts object containing singlets.
#' @param spUnsupervised An spUnsupervised object.
#' @return A matrix containing the median value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords countsPerClass
#' @examples
#'
#' output <- countsPerClass(test_spCountsSng, test_spUnsupervised)
#'
NULL

#' @rdname countsPerClass
#' @importFrom dplyr group_by summarise right_join
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export

countsPerClass <- function(
  spCountsSng,
  spUnsupervised
){
  sampleName <- NULL
  counts <- getData(spCountsSng, "counts")
  tibble(
    class = getData(spUnsupervised, "classification"),
    sampleName = colnames(getData(spCountsSng, "counts"))
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
#' @aliases expectedInteractionFreq
#' @param spUnsupervised An spUnsupervised object.
#' @return A table with the expected frequencies..
#' @author Jason T. Serviss
#' @keywords expectedInteractionFreq
#' @examples
#'
#' expectedInteractionFreq(test_spUnsupervised)
#'
NULL

#' @rdname expectedInteractionFreq
#' @importFrom dplyr "%>%"
#' @export

expectedInteractionFreq <- function(spUnsupervised) {
  getData(spUnsupervised, "classification") %>%
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
#' @param spCountsSng An spCounts object containing singlets.
#' @param spCountsMul An spCounts object containing multiplets.
#' @return A table with the expected frequencies..
#' @author Jason T. Serviss
#' @keywords estimateTotalConnections
#' @examples
#'
#' expectedInteractionFreq(test_spUnsupervised)
#'
NULL

#' @rdname estimateTotalConnections
#' @importFrom dplyr "%>%"
#' @export

estimateTotalConnections <- function(spCountsSng, spCountsMul) {
  sampleType <- cellNumberMedian <- cellNumber <- connections <- NULL
  estimateCells(spCountsSng, spCountsMul) %>%
    filter(sampleType == "Multiplet") %>%
    mutate(cellNumber = round(cellNumberMedian)) %>%
    filter(cellNumber > 1) %>%
    mutate(connections = map_dbl(cellNumber, function(n) {
      ncol(combn(1:n, 2))
    })) %>%
    pull(connections) %>%
    sum()
}

################################################################################
#                                                                              #
#                             knnClassification                                #
#                                                                              #
################################################################################

#' runPCA
#'
#' Runs a principal component analysis using \code{\link[gmodels]{fast.prcomp}}.
#'
#' @name runPCA
#' @rdname runPCA
#' @author Jason T. Serviss
#' @param spCounts spCounts object with singlets only.
#' @param select A numeric vector indicating the indexes of genes to include.
#' @importFrom gmodels fast.prcomp
#' @export

runPCA <- function(spCounts, select) {
  counts.log <- getData(spCounts, "counts.log")
  gmodels::fast.prcomp(t(counts.log[select, ]), center = TRUE)
}

#' pcByVarPercent
#'
#' Runs a principal component analysis using \code{\link[gmodels]{fast.prcomp}}.
#'
#' @name pcByVarPercent
#' @rdname pcByVarPercent
#' @author Jason T. Serviss
#' @param pca Output from \code{\link[gmodels]{fast.prcomp}} or
#' \code{\link[stats]{prcomp}}.
#' @param cutoff numeric; Select PCs containing variance greater or equal to
#' this cutoff.
#' @return The maximum PC that contains variance specified by the cutoff
#' argument.
#' @export

pcByVarPercent <- function(pca, cutoff) {
  var <- pca$sdev^2
  var.percent <- var / sum(var) * 100
  which(var.percent <= cutoff)[1] - 1
}

#' constructGraph
#'
#' Constructs a graph from a data frame or tibble.
#'
#' @name constructGraph
#' @rdname constructGraph
#' @author Jason T. Serviss
#' @param data Tibble; Expected columns are "from" and "to", indicating the
#' edges and "dist" indicating the distance between the samples.
#' @param distCut numeric; A Euclidean distance at which NN are removed.
#' @importFrom tidygraph as_tbl_graph activate .N
#' @importFrom igraph graph_from_data_frame
#' @importFrom dplyr mutate "%>%"
#' @export

constructGraph <- function(data, distCut) {
  edges <- from <- to <- nodes <- NULL

  g <- igraph::graph_from_data_frame(data, directed = FALSE) %>%
    #The function below simplifies the edges such that if node1 -> node2 and node2 -> node1, only one edge represents the connection.
    #Although, it is worth while to consider what happens if node1 -> node2 but not node2 -> node1. In a mutual kNN graph both node1 -> node2 and node2 -> node1 need to be true for an edge to be drawn between them.
    #(http://www.tml.cs.uni-tuebingen.de/team/luxburg/publications/MaiHeiLux07.pdf).
    #Consider this, if we observe the distances between NN of a small population, with a k high enough to permit NN to another class,
    #once the neighbors start to be of another class we should see a change in the distribution of the distance metric.
    igraph::simplify(edge.attr.comb = list(
      mutual = "first",
      rank = "concat",
      dist = function(d) {
        if(length(unique(d)) == 1) {unique(d)} else {d}
      }
    )) %>%
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(distOK = if_else(dist < distCut, TRUE, FALSE)) %>%
    dplyr::mutate(
      from.name = tidygraph::.N()$name[from],
      to.name = tidygraph::.N()$name[to]
    ) %>%
    tidygraph::activate(nodes)
}

#' getKNN
#'
#' Calculates Euclidean distance in PCA space.
#'
#' @name getKNN
#' @rdname getKNN
#' @author Jason T. Serviss
#' @param pca Matrix; Principal components as columns, samples as rows.
#' @param pcs Integer; Length 1 vector indicating the max principal component to
#' retain.
#' @param k Integer; Number of nearest neighbors.
#' @param mask Integer, Principal components to mask.
#' @importFrom RANN nn2
#' @importFrom tidyr gather
#' @importFrom dplyr mutate select distinct "%>%"
#' @importFrom stringr str_replace
#' @importFrom purrr map2_chr
#' @export

getKNN <- function(pca, pcs, k, mask) {
  from <- to <- x <- y <- NULL
  k2 <- k + 1
  p <- pca$x
  keep.pcs <- 1:pcs
  if(!is.null(mask)) keep.pcs <- keep.pcs[-mask]
  near_data <- RANN::nn2(p[, keep.pcs], k = k2, searchtype = 'standard', eps = 0)
  index <- near_data$nn.idx
  dists <- near_data$nn.dists
  rownames(index) <- rownames(p)
  rownames(dists) <- rownames(p)

  ds <- matrix_to_tibble(dists[, -1], "sample") %>%
    gather(rank, dist, -sample) %>%
    mutate(rank = as.numeric(str_replace(rank, "^.(.*)", "\\1")))

  matrix_to_tibble(index, drop = TRUE) %>%
    setNames(c("from", 1:k)) %>%
    tidyr::gather(rank, index, -from) %>%
    mutate(rank = as.numeric(rank)) %>%
    dplyr::mutate(to = rownames(p)[index]) %>%
    dplyr::mutate(from = rownames(p)[from]) %>%
    dplyr::select(from, to, rank) %>%
    inner_join(ds, by = c("from" = "sample", "rank" = "rank")) %>%
    #make mutual kNN
    unite(x, from, to, sep = "-", remove = FALSE) %>%
    unite(y, to, from, sep = "-", remove = FALSE) %>%
    #filter(x %in% y) %>%
    mutate(mutual = if_else(x %in% y, TRUE, FALSE)) %>%
    select(-x, -y)
}

#' jaccardSimilarity
#'
#' Calculates the Jaccard index pairwise for graph nodes..
#'
#' @name jaccardSimilarity
#' @rdname jaccardSimilarity
#' @author Jason T. Serviss
#' @param g tbl_graph; A tidygraph graph.
#' @param prune numeric; Jaccard index less than prune are marked for non inclusion
#'  in community detection.
#' @param seed integer; a seed to set before running the community detection.
#' @importFrom tidyr gather
#' @importFrom dplyr mutate "%>%"
#' @importFrom stats dist
#' @export

jaccardSimilarity <- function(g, prune = 0, seed = 7823) {
  name <- to <- from <- NULL
  set.seed(seed)
  j <- igraph::similarity(g, method = "jaccard", mode = "all")
  g <- tidygraph::activate(g, nodes)
  rownames(j) <- pull(g, name)
  colnames(j) <- pull(g, name)
  ji <- j %>%
    matrix_to_tibble("from") %>%
    gather(to, ji, -from)

  g %>%
    tidygraph::activate(edges) %>%
    inner_join(ji, by = c("from.name" = "from", "to.name" = "to")) %>%
    mutate(pruneOK = if_else(ji > prune, TRUE, FALSE)) %>%
    tidygraph::activate(nodes)
}

#' classify_louvain
#'
#' Runs the louvain community detection algorithm on the input graph.
#'
#' @name classify_louvain
#' @rdname classify_louvain
#' @author Jason T. Serviss
#' @param g tbl_graph; A tidygraph graph.
#' @param weights.col Bare word; indicating the column name of the column
#' including the weights.
#' @importFrom tidygraph group_louvain
#' @importFrom readr parse_factor
#' @importFrom dplyr mutate enquo
#' @importFrom rlang "!!"
#' @export

classify_louvain <- function(g, weights.col) {
  pruneOK <- distOK <- mutual <- louvain <- NULL
  g %>%
    tidygraph::activate(edges) %>%
    filter(pruneOK & distOK & mutual) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(
      louvain = tidygraph::group_louvain(weight = !! enquo(weights.col)),
      louvain = as.character(louvain)
    ) %>%
    tidygraph::activate(edges) %>%
    tidygraph::graph_join(g %>% tidygraph::activate(edges) %>% filter(!pruneOK | !distOK)) %>%
    tidygraph::activate(nodes)
}

#' removeOutliers
#'
#' Calculates Euclidean distance in PCA space.
#'
#' @name removeOutliers
#' @rdname removeOutliers
#' @author Jason T. Serviss
#' @param g tbl_graph; A tidygraph graph.
#' @param classificationMethod Bare word; Name of the column containing the
#' classification results
#' @param threshold Integer; Where classes containing a number of cells equal to
#'  or less than the threshold will be grouped with the most similar class.
#' @importFrom dplyr enquo count filter inner_join quo_name mutate if_else slice select pull
#' @importFrom rlang "!!" sym
#' @importFrom purrr map2_chr map2_lgl
#' @importFrom tidygraph activate
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @export

removeOutliers <- function(g, classificationMethod, threshold = 1) {
  name <- NULL
  g <- tidygraph::activate(g, nodes)
  tabl <- g %>% as_tibble() %>% count(!! rlang::sym(classificationMethod))

  ok <- filter(tabl, n > threshold) %>% inner_join(as_tibble(g), by = classificationMethod)

  undefined <- tabl %>%
    filter(n <= threshold) %>%
    inner_join(as_tibble(g), by = classificationMethod) %>%
    pull(name)

  dplyr::mutate(g, louvain = map2_chr(name, !! rlang::sym(classificationMethod), function(n, c) {
    if_else(!n %in% undefined, c, "undefined")
  }))
}
