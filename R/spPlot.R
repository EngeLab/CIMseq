
#'@include All-classes.R
NULL

#' spPlot
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name spPlot
#' @rdname spPlot
#' @aliases spPlot
#' @param counts Counts matrix with samples as columns and genes as rows.
#' @param type Can be "ercc", "markers", or .......
#' @param markers A character vector with 2 markers to plot.
#' @param ... additional arguments to pass on
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords spPlot
#' @examples
#'
#' #use demo data
#' data(Doublet_project_data)
#'
#' #run function
#' spPlot(x, type="ercc")
#'
NULL

#' @rdname spPlot
#' @export

setGeneric("spPlot", function(x, ...
){ standardGeneric("spPlot") })

#' @rdname spPlot
#' @export
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_economist

setMethod("spPlot", "spCounts", function(
    x,
    type,
    markers = NULL,
    ...
){
    if( type == "ercc" ) {
        p <- .erccPlot(x)
        p
        return(p)
    }
    if( type == "markers" ) {
        p <- .markersPlot(x, markers)
        p
        return(p)
    }
})

.erccPlot <- function(x) {
    counts <- getData(x, "counts")
    counts.ercc <- getData(x, "counts.ercc")
    sampleType <- getData(x, "sampleType")
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    d <- data.frame(sampleType = sampleType, frac.ercc=frac.ercc)
    
    p <- ggplot(d, aes(x=sampleType, y=frac.ercc))+
    geom_jitter()+
    labs(
        x="Sample type",
        y="Fraction of ERCC",
        title="Fraction ERCC in singlets/doublets"
    )+
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

.markersPlot <- function(x, markers) {
    counts.log <- getData(x, "counts.log")
    groups <- getData(x, "sampleType")
    d <- data.frame(
        sampleType = groups,
        marker1 = counts.log[markers[1], ],
        marker2 = counts.log[markers[2], ]
    )
    
    colors <- .setColors()[1:2]
    p <- ggplot(d, aes(x=marker1, y=marker2, colour=sampleType))+
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


#' @rdname spPlot
#' @export
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom reshape2 melt

setMethod("spPlot", "spUnsupervised", function(
    x,
    type,
    markers = NULL,
    ...
){
    if( type == "clusters" ) {
        p <- .spClustersPlot(x)
        p
        return(p)
    }
    if( type == "markers" ) {
        p <- .spMarkersPlot(x, markers)
        p
        return(p)
    }
})

.spClustersPlot <- function(x) {
    
    tsne <- getData(x, "tsne")
    mclust <- getData(x, "mclust")
    classification <- mclust$classification
    
    d <- cbind(as.data.frame(tsne[ ,1:2]), classification=classification)
    colors <- .setColors()
    
    p <- ggplot(d, aes(x=V1, y=V2, colour=factor(classification)))+
    geom_point(size=3, alpha=0.75)+
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
}

.spMarkersPlot <- function(x, markers) {
    tsne <- as.data.frame(getData(x, "tsne"))
    counts.log <- getData(x, "counts.log")
    sng <- counts.log[ , colnames(getData(x ,"dist"))]
    
    markExpress <- t(sng[rownames(sng) %in% markers, ])
    markExpress <- apply(markExpress, 2, function(x) {
        (x-min(x))/(max(x)-min(x))
    })
    
    df <- cbind(tsne, markExpress)
    df$sample <- rownames(df)
    m <- melt(df, id.vars=c("V1", "V2", "sample"))
    ggplot(m, aes(x=V1, y=V2, colour=value))+
        geom_point()+
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
        guides(colour=guide_colourbar(barwidth=20))
        
}


#.spMarkersPlot <- function(x, markers) {
#    markers <- c("EPCAM", "THY1", "FLT1", "INS", "GCG", "NEUROG3", "PROM1", "SST")
#
#    tsne <- getData(x, "tsne")
#    counts.log <- getData(expData, "counts.log")
#    sampleType <- getData(expData, "sampleType")
#
#    colors <- .setColors()
#    targets <- colors[1:length(markers)]
#    values <- counts.log[markers, ]
#
#    cols <- .col.from.targets(targets, values)
#
#    cex=1
#    rim.modifier=0.85
#    plot(tsne, col="black", pch=19, cex=cex)
#    points(tsne, col=cols, pch=19, cex=cex*rim.modifier)
#    legend("bottomright", legend=markers, fill=colors)
#}

.col.from.targets <- function(targets, values) {
    
    v <- t(apply(values, 1, function(x) {(x-min(x))/(max(x)-min(x))}))
    fractions <- apply(v, 2, function(x) {x/sum(x)})
    fractions[is.nan(fractions)] <- 1.0/(dim(fractions)[1])
    
    targets.rgb <- col2rgb(targets)
    res <- vector("character", length=length(targets))
    
    for(i in 1:(dim(values)[2])) {
        mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(x) {(255-x) * v[,i]}))
        mytarget.rgb <- rowSums(t(apply(mytarget.rgb, 1, function(x) {x * fractions[,i]})))
        res[i] <- rgb(red=mytarget.rgb['red']/256, green=mytarget.rgb['green']/256, blue=mytarget.rgb['blue']/256)
    }
    return(res)
}



.setColors <- function() {
    colors <- c('#0075DC', '#005C31', '#4C005C', '#2BCE48', '#FFCC99', '#808080', '#94FFB5', '#8F7C00', '#9DCC00', '#C20088', '#003380', '#FFA405', '#FFA8BB', '#426600', '#FF0010', '#5EF1F2', '#00998F', '#E0FF66', '#740AFF', '#990000', '#FFFF80', '#FFFF00', '#FF5005', '#993D59', '#00FFEF', '#FF6440', '#CC1D14', '#40FF6C', '#893D99', '#4800FF', '#FFDD40', '#CC9914', '#993D3D', '#CCCA14', '#3D996E', '#993D3D', '#191919')
    return(colors)
}

#' @rdname spPlot
#' @export
#' @import ggraph
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom igraph graph_from_data_frame V vertex
#' @importFrom ggforce theme_no_axes

setMethod("spPlot", "spSwarm", function(
    x,
    ...
){
    #get data
    spUnsupervised <- getData(x, "spUnsupervised")
    tsneMeans <- getData(spUnsupervised, "tsneMeans")
    colnames(tsneMeans) <- c("name", "x", "y")
    tsneMeans$name <- as.character(tsneMeans$name)

    tsne <- getData(spUnsupervised, "tsne")
    classification <- getData(spUnsupervised, "mclust")$classification
    d <- cbind(as.data.frame(tsne[ ,1:2]), classification=classification, stringsAsFactors=FALSE)
    colnames(d) <- c("x", "y", "name")
    
    codedSwarm <- getData(x, "codedSwarm")
    
    #process
    graph <- .networkDF(codedSwarm)
    
    #calculate the edge weights (frequency)
    graph <- .calculateWeights(graph)
    
    
    #convert to igraph
    graphDF <- graph_from_data_frame(graph)
    vertexToAdd <- tsneMeans$name[!tsneMeans$name %in% V(graphDF)$name]
    
    if(length(vertexToAdd) > 0) {
        for(i in 1:length(vertexToAdd)) {
            graphDF <- graphDF + vertex(as.character(vertexToAdd[i]))
        }
    }
    
    #setup colors
    colors <- .gg_color_hue(length(unique(tsneMeans$name)))
    names(colors) <- tsneMeans$name
    
    #setup layout
    layout <- createLayout(graphDF, 'manual', node.positions = tsneMeans)
    layout$x <- tsneMeans$x[sapply(layout$name, function(x) which(tsneMeans$name == x))]
    layout$y <- tsneMeans$y[sapply(layout$name, function(x) which(tsneMeans$name == x))]
    
    #plot
    plot <- ggraph(graph = graphDF, layout = 'manual', node.positions = layout)+
        geom_edge_link(edge_colour="black", aes(edge_alpha=weight))+
        geom_edge_loop(
            edge_colour="black",
            aes(
                angle=90,
                direction=270,
                strength=50,
                edge_alpha=weight
            )
        )+
        scale_colour_manual(values=colors)+
        geom_node_point(data=d, aes(colour=name), alpha=0.25)+
        geom_node_point(aes(colour=name), size=5)
        
    plot
    
    return(plot)
})

.gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

.networkDF <- function(x) {
    names <- colnames(x)[c(-1,-2)]
    for(o in 1:nrow(x)) {
        ind <- which(x[o, c(-1,-2)] != 0)
        
        if(length(ind) == 1) {
            combs <-  data.frame(V1=names[ind], V2=names[ind])
        } else {
            combs <- as.data.frame(t(combn(names[ind],2)), stringsAsFactors=FALSE)
        }
        
        if( o == 1 ) {
            connections <- combs
        } else {
            connections <- rbind(connections, combs)
        }
    }
    
    colnames(connections) <- c("from", "to")
    return(connections)
}

.calculateWeights <- function(graph) {
    t <- .squareTable(graph$from, graph$to)
    col <- colnames(t)
    row <- rownames(t)
    
    for(u in 1:nrow(graph)) {
        x <- as.character(graph[u, "from"])
        y <- as.character(graph[u, "to"])
        
        if(x %in% col & y %in% row) {
            graph[u, "weight"] <- t[y, x]
        } else if(y %in% col & x %in% row) {
            graph[u, "weight"] <- t[x, y]
        } else {
            graph[u, "weight"] <- 0
        }
    }
    
    graph$weight[graph$from == graph$to] <- 1
    graph$weight <- factor(graph$weight)
    return(unique(graph))
}

.squareTable <- function(x,y) {
    x <- factor(x)
    y <- factor(y)
    
    commonLevels <- sort(unique(c(levels(x), levels(y))))
    
    x <- factor(x, levels = commonLevels)
    y <- factor(y, levels = commonLevels)
    
    tbl <- table(x,y)
    t <- abs(tbl+t(tbl))
    return(t)
}



