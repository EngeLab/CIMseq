#' plotCountsData
#'
#' Assembles all data for plotCounts plots.
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
#' plotCountsData(cObjSng, cObjMul)
NULL

#' @rdname plotCountsData
#' @export

setGeneric("plotCountsData", function(
  spCountsSng,
  ...
){
  standardGeneric("plotCountsData")
})

#' @rdname plotCountsData
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%" rename mutate select full_join
#' @importFrom readr parse_factor

setMethod("plotCountsData", "spCounts", function(
  spCountsSng,
  spCountsMul,
  markers = NULL,
  ...
){
  s <- getData(spCountsSng, "counts.log")
  m <- getData(spCountsMul, "counts.log")
  
  .processMarkers(cbind(s, m), markers) %>%
  full_join(
    estimateCells(spCountsSng, spCountsMul),
    by = c("Sample" = "sampleName")
  ) %>%
  mutate(`Sample type` = parse_factor(
    sampleType,
    levels = c("Singlet", "Multiplet")
  )) %>%
  rename(`Cell number` = cellNumberMedian) %>%
  select(Sample, `Sample type`, frac.ercc, `Cell number`, 2:3)
})

#get and process data for markers plot
.processMarkers <- function(counts.log, markers) {
  
  if(is.null(markers)) {
    return(tibble(Sample = colnames(counts.log)))
  }
  
  #check that specified markers exist in data
  if(!all(markers %in% rownames(counts.log))) {
    notFound <- markers[!markers %in% rownames(counts.log)]
    notFound <- paste(notFound, collapse = ", ")
    message <- "These markers were not found in the dataset:"
    stop(paste(message, notFound))
  }
  
  #normalize the marker expression
  markExpress <- t(counts.log[rownames(counts.log) %in% markers, ])
  
  if(length(markers) == 1) {
    markExpressNorm <- matrix(normalizeVec(markExpress), ncol = 1)
  } else {
    markExpressNorm <- apply(markExpress, 2, normalizeVec)
  }
  
  #tidy markers
  markExpressNorm %>%
  matrix_to_tibble(., rowname = "Sample")
}

#' plotUnsupervisedData
#'
#' Assembles all data for plotUnsupervised plots.
#'
#' @name plotUnsupervisedData
#' @rdname plotUnsupervisedData
#' @aliases plotUnsupervisedData
#' @param spUnsupervised spUnsupervised; An spUnsupervised object.
#' @param spCountsSng spCounts; An spCounts object containing singlets.
#' @param ... additional arguments to pass on.
#' @return A tibble with columns:
#' @author Jason T. Serviss
#' @keywords plotUnsupervisedData
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' uObj <- testUns
#'
#' plotUnsupervisedData(sObjSng, uObj)
NULL

#' @rdname plotUnsupervisedData
#' @export

setGeneric("plotUnsupervisedData", function(
    spUnsupervised,
    ...
){
    standardGeneric("plotUnsupervisedData")
})

#' @rdname plotUnsupervisedData
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%" rename group_by ungroup mutate arrange summarize select full_join
#' @importFrom tibble tibble rownames_to_column as_tibble add_column
#' @importFrom tidyr gather unnest spread
#' @importFrom purrr pmap
#' @importFrom grDevices col2rgb
#' @importFrom readr parse_factor

setMethod("plotUnsupervisedData", "spUnsupervised", function(
  spUnsupervised,
  spCountsSng,
  markers = NULL,
  pal = NULL,
  ...
){
  tidyUnsupervised(spUnsupervised) %>%
  #add colors
  full_join(
    .col.from.targets(pal, getData(spCountsSng, "counts"), markers),
    by = "Sample"
  ) %>%
  #add marker data
  full_join(
    .processMarkers(getData(spCountsSng, "counts.log"), markers),
    by = "Sample"
  )
})

#pal: a colour palette
#targets: gene names
#values: gene counts
.col.from.targets <- function(
  pal,
  counts,
  markers,
  ...
){
  if(is.null(markers) | is.null(pal)) {
    samples <- colnames(counts)
    return(tibble(Sample = samples))
  }
  
  markers <- sort(markers)
  pal <- pal[1:length(markers)]
  
  counts[rownames(counts) %in% markers, ] %>%
  matrix_to_tibble(., "geneName") %>%
  gather(Sample, count, -geneName) %>%
  #normalize
  group_by(geneName) %>%
  mutate(normalized = normalizeVec(count)) %>%
  ungroup() %>%
  #calculate fraction
  group_by(Sample) %>%
  mutate(fraction = normalized / sum(normalized)) %>%
  mutate(fraction = if_else(is.nan(fraction), 1 / n(), fraction)) %>%
  #setup initial hex colors
  group_by(Sample) %>%
  arrange(geneName) %>%
  mutate(colint = pal) %>%
  ungroup() %>%
  #convert to rgb and calculate new colors
  mutate(rgb = pmap(list(colint, normalized, fraction), function(x, y, z) {
    (255 - ((255 - col2rgb(x)) * y)) * z
  })) %>%
  unnest() %>%
  add_column(col = rep(c("r", "g", "b"), nrow(.) / 3)) %>%
  group_by(Sample, col) %>%
  summarize(sumRGB = sum(rgb) / 256) %>%
  ungroup() %>%
  spread(col, sumRGB) %>%
  #convert back to hex
  mutate(Colour = pmap_chr(list(r, g, b), function(x, y, z) {
    rgb(red = x, green = y, blue = z)
  })) %>%
  select(-(b:r)) %>%
  #fix factor levels so ggplot legend will cooperate
  #https://community.rstudio.com/t/drop-false-with-scale-fill-identity/5163/2
  mutate(Colour = parse_factor(
    Colour,
    levels = unique(c(Colour, pal[!pal %in% Colour]))
  ))
}

