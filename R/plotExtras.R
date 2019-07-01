#'@include All-classes.R
NULL

#' plotData
#'
#' Returns the data used to build plots for spCounts, spUnsupervised, and
#' spSwarm objects.
#'
#' @name plotData
#' @rdname plotData
#' @param plot A ggplot object of class "gg" "ggplot".
#' @param ... additional arguments to pass on.
#' @return A tibble containing the plot data.
#' @author Jason T. Serviss
#' @examples
#'
#' p <- plotCountsERCC(CIMseqSinglets_test, CIMseqMultiplets_test)
#' plotData(p)
#'
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
#' @param singlets CIMseqSinglets; An CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; An CIMseqMultiplets object.
#' @keywords convertToERCC
#' @export
#' @importFrom dplyr select filter pull "%>%"
#' @importFrom stats median

convertToERCC <- function(ercc, singlets, multiplets) {
  estimateCells(singlets, multiplets, warning = FALSE) %>%
  select(.data$sampleType, .data$frac.ercc) %>%
  filter(.data$sampleType == "Singlet") %>%
  pull(.data$frac.ercc) %>%
  median(., na.rm = TRUE) %>%
  multiply_by(100) %>%
  divide_by(ercc)
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
  rs <- rowSums(counts[markers, ])
  if(any(rs == 0)) {
    g <- paste(markers[which(rs == 0)], collapse = ", ")
    mess <- paste0("The following genes have 0 counts for all samples ", g)
    stop(mess)
  }
  
  markers <- sort(markers)
  pal <- pal[1:length(markers)]
  
  .f1 <- function(x, y, z) {
   
  }
  
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
    function(x, y, z)  (255 - ((255 - col2rgb(x)) * y)) * z
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
    function(x, y, z) rgb(red = x, green = y, blue = z)
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

#' longFormConnections
#'
#' Helper function for plotSwarmCircos.
#'
#' @name longFormConnections
#' @rdname longFormConnections
#' @author Jason T. Serviss
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param theoretical.max integer; See \code{\link{estimateCells}}.
NULL

#' @rdname longFormConnections
#' @importFrom dplyr "%>%" mutate select group_by ungroup
#' @importFrom purrr pmap_chr
#' @importFrom readr parse_factor
#' @importFrom tidyr unite gather

longFormConnections <- function(
  swarm, singlets, multiplets, theoretical.max = NULL
){
  from <- to <- tmp <- connectionID <- super <- direction <- connectionName <- NULL
  fractions <- getData(swarm, "fractions")
  getEdgesForMultiplet(
    swarm, singlets, multiplets, rownames(fractions), 
    theoretical.max = theoretical.max
  ) %>%
    #add connectionID
    mutate(tmp = pmap_chr(list(sample, from, to), function(x, y, z) {
      paste(x, sort(c(y, z)), collapse = "-")
    })) %>%
    mutate(sub = as.numeric(parse_factor(tmp, levels = unique(tmp)))) %>%
    select(-tmp) %>%
    mutate(super = 1:nrow(.)) %>%
    unite(connectionID, super, sub, sep = ".") %>%
    #reformat
    gather(direction, class, -sample, -connectionID) %>%
    select(-direction) %>%
    #add classes
    group_by(connectionID) %>%
    mutate(from = class[1], to = class[2]) %>%
    ungroup() %>%
    unite(connectionName, from, to, remove = FALSE, sep = "--") %>%
    #add fractions
    mutate(frac = map2_dbl(sample, class, ~fractions[.x, .y]))
}

#' draw_legend
#'
#' Helper function for plotSwarmCircos.
#'
#' @name draw_legend
#' @rdname draw_legend
#' @author Jason T. Serviss
#' @param l A grob containing a legend.
NULL

#' @rdname draw_legend
#' @importFrom gridBase baseViewports
#' @importFrom graphics frame

draw_legend <- function(l){
  ##https://stackoverflow.com/questions/25192838/arrange-base-plots-and-grid-tables-on-the-same-page
  frame()
  vps <- baseViewports()
  grid::pushViewport(vps$inner, vps$figure, vps$plot)
  grid::grid.draw(l)
  grid::popViewport(3)
}

#' g_legend
#'
#' Helper function for plotSwarmCircos.
#'
#' @name g_legend
#' @rdname g_legend
#' @author Jason T. Serviss
#' @param a.gplot A ggplot.
#' @keywords g_legend
NULL

#' @rdname g_legend
#' @importFrom ggplot2 ggplot_gtable ggplot_build

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

