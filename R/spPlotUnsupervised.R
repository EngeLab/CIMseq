#'@include All-classes.R
NULL

#' plotUnsupervisedData
#'
#' Assembles all data for plotUnsupervised plots.
#'
#' @name plotUnsupervisedData
#' @rdname plotUnsupervisedData
#' @aliases plotUnsupervisedData
#' @param spUnsupervised spUnsupervised; An spUnsupervised object.
#' @param spCountsSng spCounts; An spCounts object containing singlets.
#' @param markers character; A vector with the markers to plot. Must be present
#'  in rownames of counts.
#' @param pal character; The colour palette to use.
#' @param ... additional arguments to pass on.
#' @return A tibble with columns:
#' @author Jason T. Serviss
#' @keywords plotUnsupervisedData
#' @examples
#'
#' plotUnsupervisedData(test_spUnsupervised, test_spCountsSng)
#'
NULL

#' @rdname plotUnsupervisedData
#' @importFrom dplyr "%>%" full_join

plotUnsupervisedData <-  function(
  spUnsupervised,
  spCountsSng,
  markers = NULL,
  pal = NULL,
  ...
){
  tidyUnsupervised(spUnsupervised) %>%
  #add colors
  full_join(
    coloursFromTargets(pal, getData(spCountsSng, "counts.cpm"), markers),
    by = "Sample"
  ) %>%
  #add marker data
  full_join(
    processMarkers(getData(spCountsSng, "counts.log"), markers),
    by = "Sample"
  )
}

#' plotUnsupervisedClass
#'
#' Plot method for spUnsupervised objects to highlight dimensionality reduction
#' and classification results.
#'
#' @name plotUnsupervisedClass
#' @rdname plotUnsupervisedClass
#' @aliases plotUnsupervisedClass
#' @param spUnsupervised spUnsupervised; An spUnsupervised object.
#' @param spCountsSng spCounts; An spCounts object containing singlets.
#' @param ... additional arguments to pass on.
#' @return The ggplot2 object containing the plot. See examples or the plotting
#'  vignette for further help.
#' @author Jason T. Serviss
#' @keywords plotUnsupervisedClass
#' @examples
#'
#' p <- plotUnsupervisedClass(test_spUnsupervised, test_spCountsSng)
#'
NULL

#' @rdname plotUnsupervisedClass
#' @export

setGeneric("plotUnsupervisedClass", function(
  spUnsupervised,
  spCountsSng,
  ...
){
    standardGeneric("plotUnsupervisedClass")
})

#' @rdname plotUnsupervisedClass
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%"
#' @importFrom ggthemes theme_few

#Note to self: This ideally would also work to some degree if users use their own
# methods to perform DR and classification. POtentially it is sufficient in its
# present form but this might need to be considered further...
setMethod("plotUnsupervisedClass", c("spUnsupervised", "spCounts"), function(
  spUnsupervised,
  spCountsSng,
  ...
){
    plotUnsupervisedData(spUnsupervised, spCountsSng) %>%
    ggplot(aes_string(x = '`t-SNE dim 1`', y = '`t-SNE dim 2`')) +
    geom_point(aes_string(colour = 'Classification'), alpha = 0.7) +
    scale_colour_manual(values = col40()) +
    theme_few() +
    theme(legend.position = "top") +
    guides(colour = guide_legend(override.aes = list(size = 3)))
})


#' plotUnsupervisedMarkers
#'
#' Plot method for spUnsupervised objects to highlight "markers", typically
#' genes thought to be discreetly expressed in one cell type, in the context of
#' the dimensionality reduction and classification results.
#'
#' @name plotUnsupervisedMarkers
#' @rdname plotUnsupervisedMarkers
#' @aliases plotUnsupervisedMarkers
#' @param spUnsupervised spUnsupervised; An spUnsupervised object.
#' @param spCountsSng spCounts; An spCounts object containing singlets.
#' @param markers character; A vector with markers to be included plot.
#' @param pal character; A palette of colors with length(pal) = length(markers).
#' @param ... additional arguments to pass on.
#' @return A ggplot2 object with plot. See examples or the plotting vignette
#' for further help.
#' @author Jason T. Serviss
#' @keywords plotUnsupervisedMarkers
#' @examples
#'
#' #Plot
#' p <- plotUnsupervisedMarkers(
#'   test_spUnsupervised, test_spCountsSng, markers = "CD74"
#' )
NULL

#' @rdname plotUnsupervisedMarkers
#' @export

setGeneric("plotUnsupervisedMarkers", function(
  spUnsupervised,
  spCountsSng,
  ...
){
    standardGeneric("plotUnsupervisedMarkers")
})

#' @rdname plotUnsupervisedMarkers
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%"
#' @importFrom ggthemes theme_few
#' @importFrom viridis scale_colour_viridis

#Note to self: This ideally would also work to some degree if users use their own
# methods to perform DR and classification. POtentially it is sufficient in its
# present form but this might need to be considered further...
setMethod("plotUnsupervisedMarkers", c("spUnsupervised", "spCounts"), function(
  spUnsupervised,
  spCountsSng,
  markers = NULL,
  pal = col40(),
  ...
){
  
  if(is.null(markers)) {
    stop("At least one marker must be provided in the markers argument.")
  }
  if(!all(markers %in% rownames(getData(spCountsSng, "counts")))) {
    rn <- rownames(getData(spCountsSng, "counts"))
    msg <- paste0(
      "The following markers were not found in the spCountsSng object: ",
      rn[!rn %in% markers]
    )
    stop(msg)
  }
  
  pal <- pal[seq_along(markers)]
  names(pal) <- markers
  pal <- pal[order(names(pal))]
  
  p <- plotUnsupervisedData(spUnsupervised, spCountsSng, markers, pal) %>%
  ggplot(aes_string(x = '`t-SNE dim 1`', y = '`t-SNE dim 2`')) +
  theme_few() +
  theme(legend.position = "top")
  
  if(length(markers) == 1) {
    p + geom_point(aes_string(colour = markers)) +
    scale_colour_viridis(option = "E")
  } else {
    p + geom_point(colour = "black", shape = 21) +
    geom_point(aes_string(colour = 'Colour'), alpha = .85) +
    scale_colour_identity(
      "Markers", labels = names(pal), breaks = pal,
      guide = "legend", drop = FALSE
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  }
})

