
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
    type,
    markers = NULL,
    ...
){
    #x should be an spCounts object with singlets
    #y should be an spCounts object with multuplets
    #check that type is valid
    if( type == "ercc" ) {
        p <- .countsErccPlot(x, y)
        p
        return(p)
    }
    if( type == "markers" ) {
        
        if(length(markers) != 2) {
            stop("Markers must be a character vector of length = 2.")
        }
        
        genes <- unique(rownames(getData(x, "counts")))
        
        if(!all(markers %in% genes)) {
            stop("The specified markers are not in the counts matrix.")
        }
        
        p <- .countsMarkersPlot(x, y, markers)
        p
        return(p)
    }
})

#plot ercc plot
.countsErccPlot <- function(
    x,
    y,
    ...
){
    #add function for ERCC fraction conversion
    convertToERCC <- function(ercc, x, y) {
        estimateCells(x, y) %>%
            select(.data$sampleType, .data$frac.ercc) %>%
            filter(.data$sampleType == "Singlet") %>%
            pull(.data$frac.ercc) %>%
            median %>%
            `*` (100) %>%
            `/` (ercc)
    }
    
    breaks <- c(
            100,       50,        10,          5,         2.5,
              1,      0.5,      0.25,      0.125,      0.0625,
        0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125
    )
    
    p <- estimateCells(x, y) %>%
    ggplot(
        aes_string(
            x = 'factor(sampleType, levels = c("Singlet", "Multiplet"))',
            y = 'cellNumberMedian'
        )
    ) +
    geom_jitter() +
    labs(
        x = "Sample type",
        y = "Cell number",
        title = "Fraction ERCC in singlets/doublets"
    ) +
    theme_few() +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.title = element_text(size=17),
        axis.text = element_text(size=13),
        plot.title = element_text(
            hjust = 0.5,
            family = "Arial",
            face = "bold",
            size = 24,
            margin = margin(b = 15)
        )
    ) +
    guides(
        colour = guide_legend(override.aes = list(size = 5))
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        sec.axis = sec_axis(
            trans = ~ convertToERCC(., x, y),
            name = "% ERCC",
            breaks = breaks
        )
    ) +
    expand_limits(y = 0.05)
    
    return(p)
}

#get and process data for markers plot
.countsMarkersPlotProcess <- function(
    x,
    y,
    markers,
    ...
){
    counts.log <- cbind(
        getData(x, "counts.log"),
        getData(y, "counts.log")
    )
    
    groups <- c(
        rep(
            "Singlet",
            ncol(getData(x, "counts.log"))
        ),
        rep(
            "Multuplet",
            ncol(getData(y, "counts.log"))
        )
    )
    
    d <- data.frame(
        sampleType = groups,
        marker1 = counts.log[markers[1], ],
        marker2 = counts.log[markers[2], ]
    )
    return(d)
}

#plot markers plot
.countsMarkersPlot <- function(
    x,
    y,
    markers,
    ...
){
    
    d <-
    
    colors <- col64()[7:9]
    
    p <- .countsMarkersPlotProcess(x, y, markers) %>%
    ggplot(aes_string(x='marker1', y='marker2', colour='sampleType'))+
    geom_point(size=5, alpha=0.7)+
    labs(
        x=paste("log2( Normalized counts: ", markers[1], " )", sep=""),
        y=paste("log2( Normalized counts: ", markers[2], " )", sep=""),
        title="Cell Identity Markers"
    )+
    scale_colour_manual(name="sampleType", values=colors)+
    theme_few()+
    theme(
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text=element_text(size=15),
        plot.title=element_text(
            hjust=0.5,
            family="Arial",
            face="bold",
            size=24,
            margin=margin(b=15)
        )
    )+
    guides(
        colour=guide_legend(override.aes=list(size=5))
    )
    
    return(p)
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
#' @param x An spUnsupervised object.
#' @param y An spCounts object containing singlets.
#' @param type Can be "clusters", "markers", or .......
#' @param markers A character vector with markers to plot.
#' @param plotUncertainty Logical indicating if uncertainty should be
#'    represented as the point size in the plots.
#' @param ... additional arguments to pass on.
#' @return The plotUnsupervised function returns an object of class spCounts.
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
    x,
    ...
){
    standardGeneric("plotUnsupervised")
})

#' @rdname plotUnsupervised
#' @export
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr pull

setMethod("plotUnsupervised", "spUnsupervised", function(
    x,
    y = NULL,
    type = "clusters",
    markers = NULL,
    plotUncertainty = TRUE,
    ...
){
    #x should be a spUnsupervised object.
    #y should be null or an spCounts object containig singlets.
    
    #check that type is valid
    if( type == "" ) {stop("The type argument was not specified.")}
    if( type == "clusters" ) {
        p <- .unsupClustersPlot(x, plotUncertainty)
        p
        return(p)
    }
    if( type == "markers" ) {
        #check that markers are valid
        if(is.null(y)) {
            stop("This plot requires a spCounts object as the second argument.")
        }
        p <- .unsupMarkersPlot(x, y, markers, plotUncertainty)
        p
        return(p)
    }
})

#get and process data for clusters plot
.unsupClusterPlotProcess <- function(
    x,
    plotUncertainty,
    ...
){
    d <- getData(x, "tsne") %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    as_tibble() %>%
    add_column(classification = getData(x, "classification"))
    
    if(plotUncertainty) {
        d$uncertainty <- getData(x, "uncertainty")
    } else {
        d$uncertainty <- 3
    }
    
    return(d)
}

#plot clusters plot
.unsupClustersPlot <- function(
    x,
    plotUncertainty,
    ...
){
    
    d <- .unsupClusterPlotProcess(x, plotUncertainty)
    colors <- col64()
    
    p <- ggplot(d, aes_string(x='V1', y='V2', colour='classification'))+
    #geom_point(size=3, alpha=0.75)+
    labs(
        x="Dim 1",
        y="Dim 2",
        title="Clusters"
    )+
    scale_colour_manual(name="sampleType", values=colors)+
    theme_few()+
    theme(
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text=element_text(size=15),
        plot.title=element_text(
            hjust=0.5,
            family="Arial",
            face="bold",
            size=24,
            margin=margin(b=15)
        )
    )+
    guides(
        colour=guide_legend(override.aes=list(size=5))
    )
    
    if(plotUncertainty) {
        p <- p+geom_point(aes_string(size='uncertainty'), alpha=0.75)
    } else {
        p <- p+geom_point(size=3, alpha=0.75)
    }
    
    return(p)
}

#get and process data for markers plot
.unsupMarkerPlotProcess <- function(
    x,
    y,
    markers,
    plotUncertainty,
    ...
){
    tsne <- as.data.frame(getData(x, "tsne"))
    counts.log <- getData(y, "counts.log")
    
    #chaeck that specified markers exist in data
    if(!all(markers %in% rownames(counts.log))) {
      notFound <- markers[!markers %in% rownames(counts.log)]
      notFound <- paste(notFound, collapse = ", ")
      message <- "These markers were not found in the dataset:"
      stop(paste(message, notFound))
    }
    
    markExpress <- t(counts.log[rownames(counts.log) %in% markers, ])
    markExpressNorm <- apply(markExpress, 2, function(x) {
        (x - min(x)) / (max(x) - min(x))
    })
    
    markExpressNorm %>%
        as_tibble() %>%
        add_column(sample = rownames(markExpressNorm)) %>%
        add_column(
            uncertainty = if_else(
                rep(plotUncertainty, nrow(markExpressNorm)),
                getData(x, "uncertainty"),
                0.0001
            )
        ) %>%
        add_column(V1 = pull(tsne, .data$V1), V2 = pull(tsne, .data$V2)) %>%
        gather(
            "variable",
            "value",
            -.data$V1,
            -.data$V2,
            -.data$sample,
            -.data$uncertainty
        )
}

#plot markers plot
.unsupMarkersPlot <- function(
    x,
    y,
    markers,
    plotUncertainty,
    ...
){
    
    m <- .unsupMarkerPlotProcess(x, y, markers, plotUncertainty)
    
    p <- ggplot(m, aes_string(x='V1', y='V2', colour='value'))+
        facet_grid(variable~.)+
        theme_few()+
        scale_color_viridis()+
        labs(
            x="x",
            y="y",
            title="Markers",
            color="Expression"
        )+
        theme(
            legend.position="top",
            legend.title=element_blank(),
            legend.text=element_text(size=15),
            axis.title=element_text(size=17),
            axis.text=element_text(size=15),
            plot.title=element_text(
                hjust=0.5,
                family="Arial",
                face="bold",
                size=24,
                margin=margin(b=15)
            ),
            strip.text.y = element_text(size = 15)
        )+
        guides(
            colour=guide_colourbar(barwidth=20)
        )
    
    if(plotUncertainty) {
        p <- p + geom_point(aes_string(size='uncertainty'))
    } else {
       p <- p + geom_point()
    }
        return(p)
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

.col.from.targets <- function(
    targets,
    values,
    ...
){
    
    values <- values[rownames(values) %in% targets, ]
    pal <- col64()
    targets <- pal[1:length(targets)]
    targets <- targets[1:(dim(values)[1])]
    if(is.matrix(values)) {
        
        #normalize each genes values so that all genes aer on the same scale
        v <- t(apply(values, 1, function(x) {(x - min(x)) / (max(x) - min(x))}))
        
        #calculates the fraction that each gene is expressed in a sample
        #relative to all input genes expression
        fractions <- apply(v, 2, function(x) {x / sum(x)})
        
        #this fixes the Nan problem, due to division with 0, but results in
        #positive expression values for genes samples that had no expression...
        #Do these end up being white in the end?
        fractions[is.nan(fractions)] <- 1.0 / (dim(fractions)[1])
        
        #changes the color values from hex to rgb
        targets.rgb <- col2rgb(targets)
        
        #runs a for loop with one iteration per sample
        res <- vector("character", length=length(targets))
        for(i in 1:ncol(values)) {
            
            #for each gene in the current sample, multiply the red, green, and
            #blue values (corresponding to the current gene) with the normalized
            #expression value for that gene. Subtract this from 255 to keep the
            #color within the 0-255 range.
            mytarget.rgb <- 255-t(
                apply(
                    targets.rgb,
                    1,
                    function(x) {
                        (255-x) * v[,i]
                    }
                )
            )
            
            #multiply the rgb color values for this gene by the fraction by
            #which this gene is expressed compared to other target genes
            mytarget.rgb <- rowSums(
                t(
                    apply(
                        mytarget.rgb,
                        1,
                        function(x) {
                            x * fractions[,i]
                        }
                    )
                )
            )
            
            #convert the rgb colors back to hex and output
            res[i] <- rgb(
                red=mytarget.rgb['red']/256,
                green=mytarget.rgb['green']/256,
                blue=mytarget.rgb['blue']/256
            )
        }
        return(res)
    } else {
        v <- (values-min(values))/(max(values)-min(values))
        #        fractions <- apply(v, 2, function(values) {values/sum(values)})
        #        fractions[is.nan(fractions)] <- 1.0/(dim(fractions)[1])
        targets.rgb <- col2rgb(targets)
        res <- vector("character", length=length(targets))
        for(i in 1:length(values)) {
            mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(values)
            {(255-values) * v[i]}))
            res[i] <- rgb(
                red=mytarget.rgb[1]/256,
                green=mytarget.rgb[2]/256,
                blue=mytarget.rgb[3]/256
            )
        }
        return(res)
    }
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
