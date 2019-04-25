
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
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param select A numeric vector indicating the indexes of genes to include.
#' @importFrom gmodels fast.prcomp
#' @export

runPCA <- function(singlets, select) {
  counts.log <- getData(singlets, "counts.log")
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
#' @importFrom tidyr gather unite
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
