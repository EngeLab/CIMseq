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
#' @param x An spUnsupervised object.
#' @param y An spCounts object containing singlets.
#' @param ... additional arguments to pass on.
#' @return A tibble with columns:
#' @author Jason T. Serviss
#' @keywords plotCountsData
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#'
#' \dontrun{sp.scRNAseq:::plotCountsData(cObjSng, cObjMul)}
NULL

#' @rdname plotCountsData
#' @import ggplot2
#' @importFrom dplyr "%>%" rename mutate select full_join
#' @importFrom readr parse_factor

plotCountsData <- function(
  spCountsSng,
  spCountsMul,
  markers = NULL,
  ...
){
  s <- getData(spCountsSng, "counts.log")
  m <- getData(spCountsMul, "counts.log")
  
  processMarkers(cbind(s, m), markers) %>%
  full_join(
    estimateCells(spCountsSng, spCountsMul),
    by = c("Sample" = "sampleName")
  ) %>%
  mutate(sampleType = parse_factor(
    sampleType,
    levels = c("Singlet", "Multiplet")
  )) %>%
  rename(
    `Cell number` = cellNumberMedian,
    `Sample type` = sampleType,
    `Fraction ERCC` = frac.ercc
  ) %>%
  select(-cellNumberMin, -cellNumberMax)
}

#' plotCountsERCC
#'
#' Plot method for spCounts objects to display ERCC fraction and estimated cell
#' number per sample.
#'
#' @name plotCountsERCC
#' @rdname plotCountsERCC
#' @aliases plotCountsERCC
#' @param spCountsSng spCounts; An spCounts object containing singlets.
#' @param spCountsMul spCounts; An spCounts object containing multiplets.
#' @param ... additional arguments to pass on.
#' @return A ggplot2 object containing the plot. See examples or the plotting
#'  vignette for further help.
#' @author Jason T. Serviss
#' @keywords plotCountsERCC
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#'
#' #ERCC plot
#' p <- plotCountsERCC(cObjSng, cObjMul)
#'
NULL

#' @rdname plotCountsERCC
#' @export

setGeneric("plotCountsERCC", function(
  spCountsSng,
  spCountsMul,
  ...
){
  standardGeneric("plotCountsERCC")
})

#' @rdname plotCountsERCC
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%"
#' @importFrom ggthemes theme_few

setMethod("plotCountsERCC", c("spCounts", "spCounts"), function(
    spCountsSng,
    spCountsMul,
    ...
){
  plotCountsData(spCountsSng, spCountsMul) %>%
  ggplot(aes_string(x = "`Sample type`", y = "`Cell number`")) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ convertToERCC(., cObjSng, cObjMul),
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
#' @aliases plotCountsMarkers
#' @param spCountsSng spCounts; An spCounts object containing singlets.
#' @param spCountsMul spCounts; An spCounts object containing multiplets.
#' @param markers character; A vector with the 2 markers to plot.
#' @param ... additional arguments to pass on.
#' @return A ggplot2 object containing the plot. See examples or the plotting
#'  vignette for further help.
#' @author Jason T. Serviss
#' @keywords plotCountsMarkers
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#'
#' #Markers plot
#' markers <- c("a1", "a10")
#' p <- plotCountsMarkers(cObjSng, cObjMul, markers = markers)
NULL

#' @rdname plotCountsMarkers
#' @export

setGeneric("plotCountsMarkers", function(
  spCountsSng,
  spCountsMul,
  ...
){
  standardGeneric("plotCountsMarkers")
})

#' @rdname plotCountsMarkers
#' @export
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_economist

setMethod("plotCountsMarkers", c("spCounts", "spCounts"), function(
    spCountsSng,
    spCountsMul,
    markers = NULL,
    ...
){
  if((!is.null(markers)) & length(markers) != 2) {
    stop("Markers must be a character vector of length = 2.")
  }
  
  plotCountsData(spCountsSng, spCountsMul, markers) %>%
  ggplot(aes_string(markers[1], markers[2], colour = "`Sample type`")) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = c("#1c54a8", "#f63b32")) +
  theme_few() +
  theme(legend.position = "top")
})



