#'@include All-classes.R
NULL

#' plotCountsData
#'
#' Assembles all data for plotCounts plots. Not exported. End users should use
#' the \code{\link{plotData}} function.
#'
#' @name plotCountsData
#' @rdname plotCountsData
#' @aliases plotCountsData
#' @param singlets CIMseqSinglets; An CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; An CIMseqMultiplets object.
#' @param markers character; A vector with the 2 markers to plot. Must be
#'  present in rownames of counts.
#' @param ... additional arguments to pass on.
#' @return A tibble with columns:
#' @author Jason T. Serviss
#' @keywords plotCountsData
#' @examples
#' d <- plotCountsData(CIMseqSinglets_test, CIMseqMultiplets_test)
NULL

#' @rdname plotCountsData
#' @export

setGeneric("plotCountsData", function(
  singlets, multiplets, ...
){
  standardGeneric("plotCountsData")
})

#' @rdname plotCountsData
#' @import ggplot2
#' @importFrom dplyr "%>%" rename mutate select full_join
#' @importFrom rlang .data
#' @importFrom readr parse_factor

setMethod("plotCountsData", c("CIMseqSinglets", "CIMseqMultiplets"), function(
  singlets, multiplets, markers = NULL, ...
){
  
  s <- getData(singlets, "counts.log")
  m <- getData(multiplets, "counts.log")
  
  processMarkers(cbind(s, m), markers) %>%
  full_join(
    estimateCells(singlets, multiplets),
    by = c("Sample" = "sampleName")
  ) %>%
  mutate(sampleType = parse_factor(
    .data$sampleType,
    levels = c("Singlet", "Multiplet")
  )) %>%
  full_join(., tibble(
    Sample = colnames(getData(singlets, "counts")),
    Classification = getData(singlets, "classification")
  ), by = "Sample") %>%
  full_join(
    ., matrix_to_tibble(getData(singlets, "dim.red"), "Sample"), by = "Sample"
  ) %>%
  rename(
    `Cell number` = .data$cellNumberMedian,
    `Sample type` = .data$sampleType,
    `Fraction ERCC` = .data$frac.ercc,
    `dim.red dim 1` = .data$V1,
    `dim.red dim 2` = .data$V2
  ) %>%
  select(-.data$cellNumberMin, -.data$cellNumberMax)
})

#' plotCountsERCC
#'
#' Plot method for spCounts objects to display ERCC fraction and estimated cell
#' number per sample.
#'
#' @name plotCountsERCC
#' @rdname plotCountsERCC
#' @param singlets CIMseqSinglets; An CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; An CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return A ggplot2 object containing the plot. See examples or the plotting
#'  vignette for further help.
#' @author Jason T. Serviss
#' @examples
#' p <- plotCountsERCC(CIMseqSinglets_test, CIMseqMultiplets_test)
NULL

#' @rdname plotCountsERCC
#' @export

setGeneric("plotCountsERCC", function(
  singlets, multiplets, ...
){
  standardGeneric("plotCountsERCC")
})

#' @rdname plotCountsERCC
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%"
#' @importFrom ggthemes theme_few

setMethod("plotCountsERCC", c("CIMseqSinglets", "CIMseqMultiplets"), function(
  singlets, multiplets, ...
){
  plotCountsData(singlets, multiplets) %>%
    ggplot(aes_string(x = "`Sample type`", y = "`Cell number`")) +
    scale_y_continuous(
      expand = c(0, 0),
      sec.axis = sec_axis(
        trans = ~ convertToERCC(., singlets, multiplets),
        name = "% ERCC"
      )
    ) +
    geom_jitter() +
    theme_few()
})

#' plotCountsMarkers
#'
#' Plot method for spCounts objects to display "markers", typically genes
#' thought to be discreetly expressed in one cell type, in all samples.
#'
#' @name plotCountsMarkers
#' @rdname plotCountsMarkers
#' @param singlets CIMseqSinglets; An CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; An CIMseqMultiplets object.
#' @param markers character; A vector with the 2 markers to plot.
#' @param ... additional arguments to pass on.
#' @return A ggplot2 object containing the plot. See examples or the plotting
#'  vignette for further help.
#' @author Jason T. Serviss
#' @examples
#' markers <- c("ACTB", "GAPDH")
#' p <- plotCountsMarkers(CIMseqSinglets_test, CIMseqMultiplets_test, markers)
NULL

#' @rdname plotCountsMarkers
#' @export

setGeneric("plotCountsMarkers", function(
  singlets, multiplets, ...
){
  standardGeneric("plotCountsMarkers")
})

#' @rdname plotCountsMarkers
#' @export
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_economist

setMethod("plotCountsMarkers", c("CIMseqSinglets", "CIMseqMultiplets"), function(
  singlets, multiplets, markers = NULL, ...
){
  if((!is.null(markers)) & length(markers) != 2) {
    stop("Markers must be a character vector of length = 2.")
  }
  
  plotCountsData(singlets, multiplets, markers) %>%
    ggplot(aes_string(markers[1], markers[2], colour = "`Sample type`")) +
    geom_point(alpha = 0.85) +
    scale_colour_manual(values = c("#1c54a8", "#f63b32")) +
    theme_few() +
    theme(legend.position = "top")
})

#' plotUnsupervisedClass
#'
#' Plot method for spUnsupervised objects to highlight dimensionality reduction
#' and classification results.
#'
#' @name plotUnsupervisedClass
#' @rdname plotUnsupervisedClass
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return The ggplot2 object containing the plot. See examples or the plotting
#'  vignette for further help.
#' @author Jason T. Serviss
#' @examples
#'
#' p <- plotUnsupervisedClass(CIMseqSinglets_test, CIMseqMultiplets_test)
#'
NULL

#' @rdname plotUnsupervisedClass
#' @export

setGeneric("plotUnsupervisedClass", function(
  singlets, ...
){
  standardGeneric("plotUnsupervisedClass")
})

#' @rdname plotUnsupervisedClass
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%" filter
#' @importFrom ggthemes theme_few

#Note to self: This ideally would also work to some degree if users use their own
# methods to perform DR and classification. POtentially it is sufficient in its
# present form but this might need to be considered further...
setMethod("plotUnsupervisedClass", "CIMseqSinglets", function(
  singlets, multiplets, ...
){
  `Sample type` <- NULL
  plotCountsData(singlets, multiplets) %>%
    filter(`Sample type` == "Singlet") %>%
    ggplot(aes_string(x = '`dim.red dim 1`', y = '`dim.red dim 2`')) +
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
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param markers character; A vector with markers to be included plot.
#' @param pal character; A palette of colors with length(pal) = length(markers).
#' @param ... additional arguments to pass on.
#' @return A ggplot2 object with plot. See examples or the plotting vignette
#' for further help.
#' @author Jason T. Serviss
#' @examples
#'
#' #Plot
#' p <- plotUnsupervisedMarkers(
#'   CIMseqSinglets_test, CIMseqMultiplets_test, markers = "CD74"
#' )
NULL

#' @rdname plotUnsupervisedMarkers
#' @export

setGeneric("plotUnsupervisedMarkers", function(
  singlets, multiplets, ...
){
  standardGeneric("plotUnsupervisedMarkers")
})

#' @rdname plotUnsupervisedMarkers
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%" filter
#' @importFrom ggthemes theme_few
#' @importFrom viridis scale_colour_viridis

#Note to self: This ideally would also work to some degree if users use their own
# methods to perform DR and classification. POtentially it is sufficient in its
# present form but this might need to be considered further...
setMethod("plotUnsupervisedMarkers", c("CIMseqSinglets", "CIMseqMultiplets"), function(
  singlets, multiplets, markers = NULL, pal = col40(), ...
){
  `Sample type` <- NULL
  if(is.null(markers)) {
    stop("At least one marker must be provided in the markers argument.")
  }
  if(!all(markers %in% rownames(getData(singlets, "counts")))) {
    rn <- rownames(getData(singlets, "counts"))
    msg <- paste0(
      "The following markers were not found in the spCountsSng object: ",
      rn[!rn %in% markers]
    )
    stop(msg)
  }
  
  pal <- pal[seq_along(markers)]
  names(pal) <- markers
  pal <- pal[order(names(pal))]
  
  if(length(markers) == 1) {
    p <- plotCountsData(singlets, multiplets, markers, pal) %>%
      filter(`Sample type` == "Singlet") %>%
      ggplot(aes_string(x = '`dim.red dim 1`', y = '`dim.red dim 2`')) +
      theme_few() +
      theme(legend.position = "top") + 
      geom_point(aes_string(colour = markers)) +
      scale_colour_viridis(option = "E")
  } else {
    p <- plotCountsData(singlets, multiplets, markers, pal) %>%
      filter(`Sample type` == "Singlet") %>%
      full_join(
        coloursFromTargets(pal, getData(singlets, "counts.cpm"), markers),
        by = "Sample"
      ) %>%
      ggplot(aes_string(x = '`dim.red dim 1`', y = '`dim.red dim 2`')) +
      theme_few() +
      theme(legend.position = "top") +
      geom_point(colour = "black", shape = 21) +
      geom_point(aes_string(colour = 'Colour'), alpha = .85) +
      scale_colour_identity(
        "Markers", labels = names(pal), breaks = pal,
        guide = "legend", drop = FALSE
      ) +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  }
  
  p
  return(p)
})
