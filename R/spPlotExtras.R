#'@include All-classes.R
NULL

#' plotData
#'
#' Returns the data used to build plots for spCounts, spUnsupervised, and
#' spSwarm objects.
#'
#' @name plotData
#' @rdname plotData
#' @aliases plotData
#' @param plot A ggplot object of class "gg" "ggplot".
#' @param ... additional arguments to pass on.
#' @return A tibble containing the plot data.
#' @author Jason T. Serviss
#' @keywords plotData
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#'
#' #make plot
#' p <- plotCountsERCC(cObjSng, cObjMul)
#'
#' #get data
#' plotData(p)
NULL

#' @rdname plotData
#' @export

setGeneric("plotData", function(
  plot,
  ...
){
  standardGeneric("plotData")
})

#' @rdname plotData
#' @export
#' @importFrom tibble as_tibble

setMethod("plotData", "gg", function(
  plot,
  ...
){
  if(any(grepl("ggraph", class(plot[[1]])))) {
    attr(plot[[1]], "graph")
  } else {
    as_tibble(plot[[1]]) #have a look at ggplot::fortify -> broom package
  }
})

#' convertToERCC
#'
#' A function to facilitate calculation of the second axis of the plotCounts
#' type "ercc" plot.
#'
#' @name convertToERCC
#' @rdname convertToERCC
#' @author Jason T. Serviss
#' @param ercc The left axis values. Passes as ".".
#' @param spCountsSng spCounts; An spCounts object with singlets.
#' @param spCountsMul spCounts; An spCounts object with multiplets.
#' @keywords convertToERCC
#' @export
#' @importFrom dplyr select filter pull "%>%"
#' @importFrom stats median

convertToERCC <- function(ercc, spCountsSng, spCountsMul) {
  estimateCells(spCountsSng, spCountsMul, warning = FALSE) %>%
  select(.data$sampleType, .data$frac.ercc) %>%
  filter(.data$sampleType == "Singlet") %>%
  pull(.data$frac.ercc) %>%
  median(., na.rm = TRUE) %>%
  `*` (100) %>%
  `/` (ercc)
}

#' coloursFromTargets
#'
#' Diverging color palette with 40 colors.
#'
#' @name coloursFromTargets
#' @rdname coloursFromTargets
#' @author Jason T. Serviss
#' @param pal character; A colour palette with length = length(markers).
#' @param counts matrix; A matrix containing counts.
#' @param markers character; The markers to evaluate. Must be present in
#'  rownames(counts).
#' @param ... additional arguments to pass on.
#' @keywords coloursFromTargets
NULL

#' @rdname coloursFromTargets
#' @importFrom dplyr "%>%" group_by ungroup mutate arrange summarize select if_else n
#' @importFrom rlang .data
#' @importFrom tibble tibble add_column
#' @importFrom tidyr gather unnest spread
#' @importFrom purrr pmap pmap_chr
#' @importFrom grDevices col2rgb rgb
#' @importFrom readr parse_factor

coloursFromTargets <- function(
  pal,
  counts,
  markers,
  ...
){
  
  if(is.null(markers) | is.null(pal) | length(markers) == 1) {
    return(tibble('Sample' = colnames(counts)))
  }

  markers <- sort(markers)
  pal <- pal[1:length(markers)]
  
  counts[rownames(counts) %in% markers, ] %>%
  matrix_to_tibble(., 'geneName') %>%
  gather('Sample', 'count', -.data$geneName) %>%
  #normalize
  group_by(.data$geneName) %>%
  mutate('normalized' = normalizeVec(.data$count)) %>%
  ungroup() %>%
  #calculate fraction
  group_by(.data$Sample) %>%
  mutate('fraction' = .data$normalized / sum(.data$normalized)) %>%
  mutate('fraction' = if_else(is.nan(.data$fraction), 1 / n(), .data$fraction)) %>%
  #setup initial hex colors
  arrange(.data$geneName) %>%
  mutate('colint' = pal) %>%
  ungroup() %>%
  #convert to rgb and calculate new colors
  mutate('rgb' = pmap(
    list(.data$colint, .data$normalized, .data$fraction),
    function(x, y, z) {
      (255 - ((255 - col2rgb(x)) * y)) * z
    }
  )) %>%
  unnest() %>%
  add_column('col' = rep(c("r", "g", "b"), nrow(.) / 3)) %>%
  group_by(.data$Sample, .data$col) %>%
  summarize('sumRGB' = sum(.data$rgb) / 256) %>%
  ungroup() %>%
  spread('col', 'sumRGB') %>%
  #convert back to hex
  mutate('Colour' = pmap_chr(
    list(.data$r, .data$g, .data$b),
    function(x, y, z) {
      rgb(red = x, green = y, blue = z)
    }
  )) %>%
  select(-(.data$b:.data$r)) %>%
  #fix factor levels so ggplot legend will cooperate
  #https://community.rstudio.com/t/drop-false-with-scale-fill-identity/5163/2
  mutate('Colour' = parse_factor(
    .data$Colour,
    levels = unique(c(.data$Colour, pal[!pal %in% .data$Colour]))
  ))
}

#' col40
#'
#' Diverging color palette with 40 colors.
#'
#' @name col40
#' @rdname col40
#' @author Jason T. Serviss
#' @keywords col40
#' @examples
#'
#' col40()
#'
#' @export
NULL

col40 <- function() {
  c(
  "#1c54a8", "#f63b32", "#00C2A0", "#FFDBE5", "#e0c48c", "#63FFAC", "#663000",
  "#e0a81c", "#385438", "#609060", "#6A3D9A", "#548495", "#A30059", "#8FB0FF",
  "#997D87", "#4FC601", "#8ca8c4", "#3B5DFF", "#BA0900", "#DDEFFF", "#7B4F4B",
  "#A1C299", "#0AA6D8", "#00846F", "#FFB500", "#C2FFED", "#A079BF", "#C0B9B2",
  "#C2FF99", "#00489C", "#6F0062", "#EEC3FF", "#922329", "#FFF69F", "#FF8A9A",
  "#B05B6F", "#7900D7", "#BC23FF", "#9B9700", "#0089A3"
  )
}

#' processMarkers
#'
#' Helper function for plotCountsMarkers and plotUnsupervisedMarkers. Gathers
#' and returns data in an expected format.
#'
#' @name processMarkers
#' @rdname processMarkers
#' @author Jason T. Serviss
#' @param counts.log matrix; A matrix containing log2(cpm).
#' @param markers character; The markers to process. Must be present in
#'  rownames(counts.log).
#' @keywords processMarkers
NULL

#' @rdname processMarkers
#' @importFrom dplyr "%>%"
#' @importFrom tibble tibble

processMarkers <- function(counts.log, markers) {
  
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
    markExpressNorm <- matrix(
      normalizeVec(markExpress),
      ncol = 1,
      dimnames = list(colnames(counts.log), markers)
    )
  } else {
    markExpressNorm <- apply(markExpress, 2, normalizeVec)
  }
  
  #tidy markers
  matrix_to_tibble(markExpressNorm, rowname = "Sample")
}

