#'@include All-classes.R
NULL

#' plotCounts
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name plotCounts
#' @rdname plotCounts
#' @aliases plotCounts
#' @param x An spCounts object containing singlets.
#' @param y An spCounts object containing multiplets.
#' @param type Can be "ercc", "markers", or .......
#' @param markers A character vector with 2 markers to plot.
#' @param ... additional arguments to pass on.
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords plotCounts
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#'
#' #ERCC plot
#' p <- plotCounts(cObjSng, cObjMul, type = "ercc")
#'
#' #Markers plot
#' markers <- c("a1", "a10")
#' p <- plotCounts(cObjSng, cObjMul, type = "markers", markers = markers)
#'
NULL

#' @rdname plotCounts
#' @export

setGeneric("plotCounts", function(
    x,
    ...
){
    standardGeneric("plotCounts")
})

#' @rdname plotCounts
#' @export
#' @import ggplot2
#' @importFrom dplyr filter pull
#' @importFrom stringr %>%
#' @importFrom stats median
#' @importFrom ggthemes theme_few scale_colour_economist

setMethod("plotCounts", "spCounts", function(
    x,
    y,
    type = "ercc",
    markers = NULL,
    ...
){
  
  #x should be an spCounts object with singlets
  #y should be an spCounts object with multuplets
  if((!is.null(markers)) & length(markers) != 2) {
    stop("Markers must be a character vector of length = 2.")
  }
  if(!type %in% c("ercc", "markers")) {
    stop("The type argument must be ercc or markers")
  }
  
  pData <- plotCountsData(x, y, markers)
  
  if(type == "ercc") {
    pData %>%
    ggplot(aes(x = `Sample type`, y = `Cell number`))
  } else {
    pData %>%
    ggplot(aes_string(markers[1], markers[2], colour = "`Sample type`"))
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
#' @param spCountsSng spCounts; An spCounts object with multiplets.
#' @keywords convertToERCC
#'
#' @export
#' @importFrom dplyr select filter pull
#' @importFrom stats median

convertToERCC <- function(ercc, spCountsSng, spCountsMul) {
  estimateCells(spCountsSng, spCountsMul) %>%
  select(.data$sampleType, .data$frac.ercc) %>%
  filter(.data$sampleType == "Singlet") %>%
  pull(.data$frac.ercc) %>%
  median %>%
  `*` (100) %>%
  `/` (ercc)
}


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
#' #
#'
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
  #!!!!! THIS WILL ERROR IF length(markers) == 1 !!!!!!!!!!!!!!
  markExpress <- t(counts.log[rownames(counts.log) %in% markers, ])
  markExpressNorm <- apply(markExpress, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  
  #tidy markers
  markExpressNorm %>%
  matrix_to_tibble(.) %>%
  rename(Sample = rowname)
}

#' plotUnsupervised
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name plotUnsupervised
#' @rdname plotUnsupervised
#' @aliases plotUnsupervised
#' @param spUnsupervised spUnsupervised; An spUnsupervised object.
#' @param spCounts spCounts; An spCounts object containing singlets.
#' @param type Not implemented.
#' @param markers character; A vector with markers to be included plot.
#' @param pal character; A palette of colors with length(pal) = length(markers).
#' @param ... additional arguments to pass on.
#' @return The ggplot2 object with the t-SNE results plotted on the x and y
#'  axis. The plot can be modified by adding geoms, themes, etc. in the normal
#' manner with ggplot2. See examples or the plotting vignette for further help.
#' @author Jason T. Serviss
#' @keywords plotUnsupervised
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' uObj <- testUns
#'
#' #Clusters plot
#' p <- plotUnsupervised(uObj, cObjSng, type = "clusters")
#'
#' #Markers plot
#' markers <- c("a1", "a10")
#' p <- plotUnsupervised(uObj, cObjSng, type = "markers", markers = markers)
#'
NULL

#' @rdname plotUnsupervised
#' @export

setGeneric("plotUnsupervised", function(
    spUnsupervised,
    ...
){
    standardGeneric("plotUnsupervised")
})

#' @rdname plotUnsupervised
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%"

setMethod("plotUnsupervised", "spUnsupervised", function(
  spUnsupervised,
  spCounts,
  type = "clusters",
  markers = NULL,
  pal = NULL,
  ...
){
    #x should be a spUnsupervised object.
    #y should be an spCounts object containig singlets.
    p <- plotUnsupervisedData(spUnsupervised, spCounts, markers, pal) %>%
    ggplot(aes(x = `t-SNE dim 1`, y = `t-SNE dim 2`))
    
    p
    return(p)
})

#' plotUnsupervisedData
#'
#' Assembles all data for plotUnsupervised plots.
#'
#' @name plotUnsupervisedData
#' @rdname plotUnsupervisedData
#' @aliases plotUnsupervisedData
#' @param x An spUnsupervised object.
#' @param y An spCounts object containing singlets.
#' @param ... additional arguments to pass on.
#' @return A tibble with columns:
#' @author Jason T. Serviss
#' @keywords plotUnsupervisedData
#' @examples
#' #
#'
NULL

#' @rdname plotUnsupervisedData
#' @export

setGeneric("plotUnsupervisedData", function(
    x,
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
  x,
  y,
  markers = NULL,
  pal = NULL,
  ...
){
  
  #check that if markers is ! NULL that pal is also ! NULL
  
  tidyUnsupervised(x) %>%
  #add colors
  full_join(
    .col.from.targets(pal, getData(y, "counts"), markers),
    by = "Sample"
  ) %>%
  #add marker data
  full_join(
    .processMarkers(getData(y, "counts.log"), markers),
    by = "Sample"
  )
})

#pal: a colour palette
#targets: gene names
#values: gene counts
.col.from.targets <- function(
  pal,
  values,
  markers,
  ...
){
  if(is.null(markers) | is.null(pal)) {
    samples <- colnames(values)
    return(tibble(Sample = samples))
  }
  
  markers <- sort(markers)
  pal <- pal[1:length(markers)]
  
  values[rownames(values) %in% markers, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "geneName") %>%
  as_tibble() %>%
  gather(Sample, count, -geneName) %>%
  #normalize
  group_by(geneName) %>%
  mutate(normalized = (count - min(count)) / (max(count) - min(count))) %>%
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

#' plotSwarm
#'
#' Description.
#'
#' Details.
#'
#' @name plotSwarm
#' @rdname plotSwarm
#' @aliases plotSwarm
#' @param x An spSwarm object.
#' @param y An spUnsupervised object.
#' @param z An spCounts object with singlets.
#' @param w An spCounts object with multiplets.
#' @param type Specifies the plot type. Possible values are "tsne" or "igraph"
#'    for network plots, "multiplets" or "edges" for residuals plots, "edgeBar"
#'    or "pValueBar" for bar plots and "heat" for a heatmap.
#' @param loop Logical; TRUE if "self" connections should be plotted.
#' @param markers Not currently implemented.
#' @param edge.cutoff The minimum fraction to consider (?).
#' @param min.pval Minimum p-value to report.
#' @param min.num.edges Minimum number of observed edges to report a connection.
#' @param label.cutoff The number of labels to show in the residuals plot.
#' @param ... additional arguments to pass on.
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords plotSwarm
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#' uObj <- testUns
#' sObj <- testSwa
#'
#' # tsne plot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "tsne")
#'
#' #igraph plot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "igraph")
#'
#' #multiplet residuals plot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "multiplets")
#'
#' #edge residuals plot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "edges")
#'
#' #edge barplot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "edgeBar")
#'
#' #pvalue barplot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "pValueBar")
#'
#' #heatmap plot
#' p <- plotSwarm(sObj, uObj, cObjSng, cObjMul, type = "heat")
#'
NULL

#' @rdname plotSwarm
#' @export

setGeneric("plotSwarm", function(
    x,
    ...
){
    standardGeneric("plotSwarm")
})

#' @rdname plotSwarm
#' @export
#' @import ggraph
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom igraph graph_from_data_frame delete_edges set_edge_attr
#'    get.edgelist E subgraph.edges
#' @importFrom ggforce theme_no_axes
#' @importFrom utils combn
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggrepel geom_label_repel
#' @importFrom viridis scale_fill_viridis
#' @importFrom readr parse_factor
#' @importFrom dplyr pull

setMethod("plotSwarm", "spSwarm", function(
    x,
    y,
    z = NULL,
    w = NULL,
    type = NULL,
    loop = FALSE,
    markers = NULL,
    edge.cutoff = NULL,
    min.pval = NULL,
    min.num.edges = NULL,
    label.cutoff = NULL,
    ...
){
    
    #cimSng,
    #cimMul,
    #cimClass,
    #cimSwarm,
    
    #x should be an spSwarm object
    #y should be an spUnsupervised object
    #z should be an spCounts object with singlets
    #w should be an spCounts object with multiplets
    
    if(is.null(edge.cutoff)) {edge.cutoff <- 1 / ncol(getData(x, "spSwarm"))}
    if(is.null(min.num.edges)) {min.num.edges <- 0}
    if(is.null(min.pval)) {min.pval <- 1}
    if(is.null(label.cutoff)) {label.cutoff <- 1}

    types <- c(
        "tsne",
        "igraph",
        "multiplets",
        "edges",
        "edgeBar",
        "pValueBar",
        "heat"
    )
    
    if(!type %in% types) {
        stop("Correct type not specified for plot.")
    }
    
    if(type %in% c("multiplets", "edges")) {
        plot <- .residualsPlot(
            w,
            y,
            x,
            type,
            edge.cutoff,
            min.pval,
            min.num.edges,
            label.cutoff
        )
    } else if(type %in% c("tsne", "igraph")) {
        plot <- .swarmPlot(
            x,
            y,
            z,
            type,
            loop = FALSE,
            markers = NULL,
            edge.cutoff,
            min.pval,
            min.num.edges
        )
    } else if(type %in% c("edgeBar", "pValueBar")) {
        plot <- .barPlot(
            x,
            edge.cutoff,
            min.pval,
            min.num.edges,
            type
        )
    } else {
        plot <- .heatmap(
            x,
            edge.cutoff,
            min.pval,
            min.num.edges
        )
    }
    
    plot
    return(plot)
    
})

################################################################################
#                                                                              #
# Swarm Plots                                                                  #
#                                                                              #
################################################################################

.swarmPlot <- function(
    x,
    y,
    z,
    type,
    loop=FALSE,
    markers=NULL,
    edge.cutoff,
    min.pval,
    min.num.edges,
    ...
){
    #get data
    tsneMeans <- getData(y, "tsneMeans")
    d <- .swarmTsneProcess(y)
    
    if(!is.null(markers)) {
        d$color <- .col.from.targets(markers, getData(z, "counts.log"))
    }
    
    graph <- spSwarmPoisson(
        x,
        edge.cutoff = edge.cutoff,
        min.pval = 1,
        min.num.edges = 0
    )
    
    #convert to igraph
    graphDF <- graph_from_data_frame(graph)
    
    #remove edges according to min.num.edges and min.pval
    graphDF <- subgraph.edges(
        graph = graphDF,
        eids = which(E(graphDF)$weight >= min.num.edges),
        delete.vertices = FALSE
    )
    
    graphDF <- subgraph.edges(
        graph = graphDF,
        eids = which(E(graphDF)$pval < min.pval),
        delete.vertices = FALSE
    )

    #add nodes not in graph (due to min.pval argument)
    #vertexToAdd <- tsneMeans$classification[!tsneMeans[,1] %in% V(graphDF)$name]
    
    #if(length(vertexToAdd) > 0) {
    #    for(i in 1:length(vertexToAdd)) {
    #        graphDF <- graphDF + vertex(as.character(vertexToAdd[i]))
    #    }
    #}
    
    #o <- match(V(graphDF)$name, tsneMeans$classification)
    #tsneMeans <- tsneMeans[o, ]
    
    #name edges
    graphDF <- set_edge_attr(
        graphDF,
        name = "name",
        value = paste(
            get.edgelist(graphDF)[,1],
            get.edgelist(graphDF)[,2],
            sep = "->"
        )
    )
    
    #setup colors
    colors <- col64()[1:length(unique(tsneMeans$classification))]
    names(colors) <- tsneMeans$classification
    
    if(type == "tsne" & loop == TRUE) {
        plot <- .plotTsneLoop(graphDF, tsneMeans, colors, d)
    }
    
    if(type == "tsne" & loop == FALSE) {
        plot <- .plotTsne(graphDF, tsneMeans, colors, d)
    }
    
    if(type == "igraph" & loop == TRUE) {
        plot <- .plotIgraphLoop(graphDF, colors)
    }
    
    if(type == "igraph" & loop == FALSE) {
        plot <- .plotIgraph(graphDF, colors)
    }
    
    return(plot)
}

.swarmTsneProcess <- function(
    y,
    ...
){
    getData(y, "tsne") %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(var = "sample") %>%
    as_tibble() %>%
    add_column(classification = getData(y, "classification")) %>%
    rename(x = V1, y = V2) %>%
    arrange(classification)
}

.plotTsne <- function(
    graphDF,
    tsneMeans,
    colors,
    d,
    ...
){
    
    #setup layout
    layout <- create_layout(graphDF, 'manual', node.positions = tsneMeans)
    
    #plot
    plot <- ggraph(
        graph = graphDF,
        layout = 'manual',
        node.positions = layout
    )+
    geom_edge_link(
        edge_colour = "black",
        aes_string(edge_width = 'factor(weight)'),
        edge_alpha = 0.3,
        lineend = "round"
    )+
    geom_node_point(
        data = d,
        aes_string(colour = 'classification'),
        alpha = 0.3
    )+
    geom_node_point(
        data = tsneMeans,
        aes_string(colour = 'classification'),
        size = 5
    )+
    scale_colour_manual(
        name = "classification",
        values = colors
    )+
    guides(
        colour = guide_legend(title = "Classification"),
        edge_width = guide_legend(title = "Weight")
    )+
    theme_few()
        
    return(plot)
}

.plotTsneLoop <- function(
    graphDF,
    tsneMeans,
    colors,
    d,
    ...
){
    #setup layout
    layout <- create_layout(graphDF, 'manual', node.positions = tsneMeans)
    
    #plot
    plot <- ggraph(
        graph = graphDF,
        layout = 'manual',
        node.positions = layout
    )+
    geom_edge_link(
        edge_colour = "black",
        aes_string(edge_width = 'factor(weight)'),
        edge_alpha = 0.3,
        lineend = "round"
    )+
    geom_edge_loop(
        edge_colour = "black",
        edge_alpha = 0.3,
        aes_string(
            angle = 90,
            direction = 270,
            strength = 50,
            edge_width = 'factor(weight)'
        ),
        lineend = "round"
    )+
    geom_node_point(
        data = d,
        aes_string(colour = 'classification'),
        alpha = 0.3
    )+
    geom_node_point(
        data = tsneMeans,
        aes_string(colour = 'classification'),
        size = 5
    )+
    scale_colour_manual(
        name = "classification",
        values = colors
    )+
    guides(
        colour = guide_legend(title = "Classification"),
        edge_width = guide_legend(title = "Weight")
    )+
    theme_few()
    
    return(plot)
}

.plotIgraph <- function(
    graphDF,
    colors,
    ...
){
    plot <- ggraph(
        graph = graphDF,
        layout = 'igraph',
        algorithm = 'kk'
    )+
    geom_edge_link(
        edge_colour = "black",
        aes_string(
            edge_width = 'factor(weight)'
        ),
        edge_alpha = 0.3
    )+
    geom_node_point(
        aes_string(colour = 'name'),
        size = 4
    )+
    scale_colour_manual(values = colors)+
    theme_few()+
    guides(
        colour = guide_legend(title = "Classification"),
        edge_width = guide_legend(title = "Weight")
    )
    
    return(plot)
}

.plotIgraphLoop <- function(
    graphDF,
    colors,
    ...
){
    plot <- ggraph(
        graph = graphDF,
        layout = 'igraph',
        algorithm = 'kk'
    )+
    geom_edge_link(
        edge_colour = "black",
        edge_alpha = 0.3,
        aes_string(edge_width = 'factor(weight)')
    )+
    geom_edge_loop(
        edge_colour = "black",
        edge_alpha = 0.3,
        aes_string(
            angle = 90,
            direction = 270,
            strength = 1,
            edge_width = 'factor(weight)'
        )
    )+
    geom_node_point(
        aes_string(colour = 'name'),
        size = 4
    )+
    scale_colour_manual(values = colors)+
    theme_few()+
    guides(
        colour = guide_legend(title = "Classification"),
        edge_width = guide_legend(title = "Weight")
    )+
    labs(
        x = "FIX THIS"
    )
}

################################################################################
#                                                                              #
# Residual Plots                                                               #
#                                                                              #
################################################################################

.residualsPlot <- function(
    x,
    y,
    z,
    type,
    edge.cutoff,
    min.pval,
    min.num.edges,
    label.cutoff,
    ...
){
    
    #x should be a spCounts object with multiplets
    #y should be a spUnsupervised object
    #z should be a spSwarm object
    
    if(!type %in% c("multiplets", "edges")) {stop("Plot type not supported")}
    
    resid <- calcResiduals(x, y, z, edge.cutoff = edge.cutoff)
    
    switch(
        type,
        multiplets = {
            d <- .resPlotMultiplets(resid)
        },
        edges = {
            d <- .resPlotEdge(
                z,
                edge.cutoff,
                min.num.edges,
                min.pval,
                resid
            )
        }
    )
    
    #setup labels
    mults <- unique(d$multiplet)
    idx <- lapply(1:length(mults), function(j)
        which(d$multiplet == mults[j])[1:label.cutoff]
    )
    label <- d[unlist(idx),]
    
    p <- ggplot(
        d,
        aes_string(x = 'multiplet', y = 'residuals')
    ) +
    geom_violin(
        color = "grey"
    ) +
    geom_segment(
        aes_string(
            x = 'match(multiplet, levels(multiplet)) - 0.1',
            xend = 'match(multiplet, levels(multiplet)) + 0.1',
            y = 'residuals',
            yend = 'residuals'
        ),
        color = "black"
    )+
    geom_label_repel(
        data = label,
        aes_string(
            x = 'multiplet',
            y = 'residuals',
            label = 'genes'
        ),
        min.segment.length = unit(0.1, "lines")
    )+
    theme_few()+
    theme(
        axis.text.x = element_text(angle = 90)
    )+
    labs(
        x = ifelse(type == "multiplets", "Multiplet", "Edge"),
        y = "Residuals"
    )
    
    p
    return(p)
}

.resPlotMultiplets <- function(
    resid,
    ...
){
    resid %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(genes = .data$rowname) %>%
    as_tibble() %>%
    gather("multiplet", "residuals", -.data$genes) %>%
    arrange(.data$multiplet, desc(.data$residuals)) %>%
    mutate(
        multiplet = parse_factor(
            .data$multiplet,
            levels = unique(.data$multiplet)
        )
    )
}

.resPlotEdge <- function(
    spSwarm,
    edge.cutoff,
    min.num.edges,
    min.pval,
    resid
){
    
    edges <- getMultipletsForEdge(
        spSwarm,
        edge.cutoff = edge.cutoff,
        edges = spSwarmPoisson(
            spSwarm,
            edge.cutoff = edge.cutoff,
            min.num.edges = min.num.edges,
            min.pval = min.pval
        ) %>% filter(.data$weight != 0)
    ) %>%
        mutate(connection = paste(.data$from, .data$to, sep = "-"))
    
    nUM <- pull(distinct(edges, .data$connection), .data$connection)
    sapply(1:length(nUM), function(j) {
            
        muls <- filter(edges, .data$connection == nUM[j]) %>%
            pull(.data$multiplet)
            
        if(length(muls) != 1) {
            rowSums(resid[, colnames(resid) %in% muls])
        } else {
            resid[, colnames(resid) %in% muls]
        }
    }) %>%
        as.data.frame() %>%
        setNames(nUM) %>%
        rownames_to_column() %>%
        rename(genes = .data$rowname) %>%
        as_tibble() %>%
        gather("multiplet", "residuals", -.data$genes) %>%
        arrange(.data$multiplet, desc(.data$residuals)) %>%
        mutate(
            multiplet = parse_factor(.data$multiplet, unique(.data$multiplet))
        )
}

################################################################################
#                                                                              #
# Bar Plots                                                                    #
#                                                                              #
################################################################################

.barPlot <- function(
    sObj,
    edge.cutoff,
    min.pval,
    min.num.edges,
    type
){
    input <- .processFullConnections(
        sObj,
        edge.cutoff,
        min.pval,
        min.num.edges
    )
    
    if(type == "pValueBar") {
        plot <- .pvalueBar(input)
    } else {
        plot <- .edgesBar(input)
    }
    
    return(plot)
}

.processFullConnections <- function(
    sObj,
    edge.cutoff,
    min.pval,
    min.num.edges,
    ...
){
    types <- colnames(getData(sObj, "spSwarm"))
    
    d <- data.frame(
        from = sort(rep(types, length(types))),
        to   = rep(types, length(types))
    )
    
    results <- spSwarmPoisson(
        spSwarm       = sObj,
        edge.cutoff   = edge.cutoff,
        min.pval      = min.pval,
        min.num.edges = min.num.edges
    )
    
    tmp1 <- merge(d, results, by.x = c("from", "to"), by.y = c("from", "to"))
    tmp2 <- merge(d, results, by.x = c("from", "to"), by.y = c("to", "from"))
    input <- unique(rbind(tmp1, tmp2))
    
    return(input)
}

.pvalueBar <- function(
    input,
    ...
){
    
    colors <- col64()[1:length(unique(input[,1]))]
    
    p <- ggplot(
        input,
        aes_string(
            x = 'to',
            y = '-log10(pval)'
        )
    )+
    geom_bar(
        aes_string(fill = 'to'),
        stat = "identity",
        position = position_dodge(width = 1)
    )+
    facet_grid(
        from ~ to,
        scales = "free_x"
    )+
    geom_hline(
        yintercept = -log10(0.05),
        lty = 2,
        colour = "darkgrey"
    )+
    theme_few()+
    labs(
        x = "Node 1",
        y = "-log10(p-value)"
    )+
    scale_fill_manual(values = colors)+
    guides(fill = FALSE)+
    theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
    )
    
    return(p)
}

.edgesBar <- function(
    input,
    ...
){
    colors <- col64()[1:length(unique(input[,1]))]
    
    input$pStars <- ifelse(
        input$pval < 0.05 & input$pval > 0.0049,
        "\u2605",
        ifelse(
            input$pval < 0.005 & input$pval > 0.00049,
            "\u2605\u2605",
            ifelse(
                input$pval < 0.0005,
                "\u2605\u2605\u2605",
                ""
            )
        )
    )
    
    caption <- paste(
        "\u2605 = p < 0.05",
        "\u2605\u2605 = p < 0.005",
        "\u2605\u2605\u2605 = p < 0.0005.",
        sep = "; "
    )
    
    p <- ggplot(
        input,
        aes_string(
            x = 'to',
            y = 'weight'
        )
    )+
    geom_bar(
        aes_string(
            fill = 'to'
        ),
        stat = "identity",
        position = position_dodge(width = 1)
    )+
    facet_grid(
        from ~ to,
        scales = "free_x"
    )+
    geom_label(
        aes_string(
            x = 'to',
            y = 'weight + 6',
            label = 'pStars'
        ),
        fill = "white",
        label.size = 0,
        label.padding = unit(0.01, "lines"),
        position = position_dodge(width = 0.9),
        show.legend = FALSE,
        family = "Arial Unicode MS",
        size = 3
    )+
    theme_few()+
    labs(
        y = "Number of Edges",
        caption = caption
    )+
    scale_fill_manual(
        values = colors
    )+
    guides(
        fill = FALSE
    )+
    theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0),
        strip.text.y = element_text(angle = 0,  hjust = 0),
        plot.caption = element_text(
            family = "Arial Unicode MS",
            hjust  = 0
        )
    )+
    ylim(0, max(input$weight + 10))
    
    return(p)
}

################################################################################
#                                                                              #
# Heatmap                                                                      #
#                                                                              #
################################################################################

.heatmap <- function(
    sObj,
    edge.cutoff,
    min.pval,
    min.num.edges,
    ...
){
    input <- .processFullConnections(
        sObj,
        edge.cutoff,
        min.pval,
        min.num.edges
    )
    
    input$pStars <- ifelse(
        input$pval < 0.05 & input$pval > 0.0049,
        "\u2605",
        ifelse(
            input$pval < 0.005 & input$pval > 0.00049,
            "\u2605\u2605",
            ifelse(
                input$pval < 0.0005,
                "\u2605\u2605\u2605",
                ""
            )
        )
    )
    
    caption <- paste(
        "\u2605 = p < 0.05",
        "\u2605\u2605 = p < 0.005",
        "\u2605\u2605\u2605 = p < 0.0005.",
        sep = "; "
    )
    
    plot <- ggplot(
        input,
        aes_string(
            x = 'from',
            y = 'to'
        )
    )+
    geom_tile(
        aes_string(
            fill = 'weight'
        )
    )+
    theme_few()+
    labs(
        x = "From",
        y = "To",
        caption = caption
    )+
    guides(fill = guide_legend(title = "Weight"))+
    scale_fill_viridis()+
    geom_text(
        aes_string(
            label = 'pStars'
        ),
        family = "Arial Unicode MS"
    )+
    theme(
        plot.caption = element_text(
            family = "Arial Unicode MS",
            hjust = 0
        )
    )
    
    return(plot)
    
}

#' col64
#'
#' Diverging color palette.
#'
#' @name col64
#' @rdname col64
#' @author Jason T. Serviss
#' @keywords col64
#' @examples
#'
#'cols <- col64()
#'
#' @export
NULL

col64 <- function() {
    c(
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6",
    "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43",
    "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601",
    "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0",
    "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B",
    "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500",
    "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C",
    "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
    "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
    "#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700",
    "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837",
    "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625",
    "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC",
    "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"
    )
}
