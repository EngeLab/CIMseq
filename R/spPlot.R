
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
#' @param x Counts matrix with samples as columns and genes as rows.
#' @param type Can be "ercc", "markers", or .......
#' @param markers A character vector with 2 markers to plot.
#' @param layout Specify "tsne" for overlaying connections over the tSNE or igraph for graph layout.
#' @param loop Logical; TRUE if "self" connections should be turned on.
#' @param ... additional arguments to pass on
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords spPlot
#' @examples
#'
#' #use demo data
#' data(expData)
#'
#' #run function
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
    singlets,
    multuplets,
    type,
    markers = NULL,
    ...
){
    #check that type is valid
    if( type == "ercc" ) {
        p <- .countsErccPlot(singlets, multuplets)
        p
        return(p)
    }
    if( type == "markers" ) {
        p <- .countsMarkersPlot(singlets, multuplets, markers)
        p
        return(p)
    }
})

#get and process data for ercc plot
.countsErccPlotProcess <- function(sng, mul) {
    countsSng <- getData(sng, "counts")
    countsMul <- getData(mul, "counts")
    counts <- cbind(countsSng, countsMul)
    
    counts.erccSng <- getData(sng, "counts.ercc")
    counts.erccMul <- getData(mul, "counts.ercc")
    counts.ercc <- cbind(counts.erccSng, counts.erccMul)

    sampleType <- c(rep("Singlet", ncol(countsSng)), rep("Multuplet", ncol(countsMul)))
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    d <- data.frame(sampleType = sampleType, frac.ercc=frac.ercc)
    return(d)
}

#plot ercc plot
.countsErccPlot <- function(x) {
    
    d <- .countsErccPlotProcess(x)
    
    p <- ggplot(d, aes_string(x='sampleType', y='frac.ercc'))+
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

#get and process data for markers plot
.countsMarkersPlotProcess <- function(x, markers) {
    
    counts.logSng <- getData(sng, "counts.log")
    counts.logMul <- getData(mul, "counts.log")
    counts.log <- cbind(counts.logSng, counts.logMul)
    
    groups <- c(rep("Singlet", ncol(counts.logSng)), rep("Multuplet", ncol(counts.logMul)))

    d <- data.frame(
        sampleType = groups,
        marker1 = counts.log[markers[1], ],
        marker2 = counts.log[markers[2], ]
    )
    return(d)
}

#plot markers plot
.countsMarkersPlot <- function(sng, mul, markers) {
    
    d <- .countsMarkersPlotProcess(sng, mul, markers)
    
    colors <- .setColors()[1:2]
    
    p <- ggplot(d, aes_string(x='marker1', y='marker2', colour='sampleType'))+
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
#' @importFrom viridis scale_color_viridis

setMethod("spPlot", "spUnsupervised", function(
    x,
    type,
    markers = NULL,
    ...
){
    #check that type is valid
    if( type == "clusters" ) {
        p <- .unsupClustersPlot(x)
        p
        return(p)
    }
    if( type == "markers" ) {
        #check that markers are valid
        p <- .unsupMarkersPlot(x, markers)
        p
        return(p)
    }
})

#get and process data for clusters plot
.unsupClusterPlotProcess <- function(x) {
    tsne <- getData(x, "tsne")
    classification <- getData(x, "classification")
    
    d <- cbind(as.data.frame(tsne[ ,1:2]), classification=classification)
    
    return(d)
}

#plot clusters plot
.unsupClustersPlot <- function(x) {
    
    d <- .unsupClusterPlotProcess(x)
    colors <- .setColors()

    p <- ggplot(d, aes_string(x='V1', y='V2', colour='classification'))+
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
    
    return(p)
}

#get and process data for markers plot
.unsupMarkerPlotProcess <- function(x, markers) {
    tsne <- as.data.frame(getData(x, "tsne"))
    counts.log <- getData(x, "counts.log")
    sampleType <- getData(x, "sampleType")
    sng <- counts.log[ , sampleType == "Singlet"]
    
    #on the next line you will need a check that the markers exist in the data
    markExpress <- t(sng[rownames(sng) %in% markers, ])
    markExpress <- apply(markExpress, 2, function(x) {
        (x-min(x))/(max(x)-min(x))
    })
    
    markExpress <- as.data.frame(markExpress)
    markExpress$sample <- rownames(markExpress)
    rownames(markExpress) <- 1:nrow(markExpress)
    df <- cbind(tsne, markExpress)
    m <- melt(df, id.vars=c("V1", "V2", "sample"))
    
    return(m)
}

#plot markers plot
.unsupMarkersPlot <- function(x, markers) {
    
    m <- .unsupMarkerPlotProcess(x, markers)
    
    p <- ggplot(m, aes_string(x='V1', y='V2', colour='value'))+
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
        
        return(p)
}

.setColors <- function() {
    colors <- c('#0075DC', '#005C31', '#4C005C', '#2BCE48', '#FFCC99', '#808080', '#94FFB5', '#8F7C00', '#9DCC00', '#C20088', '#003380', '#FFA405', '#FFA8BB', '#426600', '#FF0010', '#5EF1F2', '#00998F', '#E0FF66', '#740AFF', '#990000', '#FFFF80', '#FFFF00', '#FF5005', '#993D59', '#00FFEF', '#FF6440', '#CC1D14', '#40FF6C', '#893D99', '#4800FF', '#FFDD40', '#CC9914', '#993D3D', '#CCCA14', '#3D996E', '#993D3D', '#191919', 'cadetblue3', 'antiquewhite4', 'aquamarine3', 'brown4', 'darkgoldenrod3')
    return(colors)
}

#' @rdname spPlot
#' @export
#' @import ggraph
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom igraph graph_from_data_frame V vertex
#' @importFrom ggforce theme_no_axes
#' @importFrom utils combn

setMethod("spPlot", "spSwarm", function(
    x,
    layout="tsne",
    loop=TRUE,
    ...
){
    #get data
    tsneMeans <- .swarmTsneMeansProcess(x)
    d <- .swarmTsneProcess(x)
    codedSwarm <- getData(x, "codedSwarm")
    
    #process
    graph <- .swarmNetworkDF(codedSwarm)
    
    #calculate the edge weights (connection number)
    graph <- .swarmCalculateWeights(graph)
    
    
    #convert to igraph
    graphDF <- graph_from_data_frame(graph)
    vertexToAdd <- tsneMeans$name[!tsneMeans$name %in% V(graphDF)$name]
    
    if(length(vertexToAdd) > 0) {
        for(i in 1:length(vertexToAdd)) {
            graphDF <- graphDF + vertex(as.character(vertexToAdd[i]))
        }
    }
    
    #setup colors
    colors <- .setColors()[1:length(unique(tsneMeans$name))]
    names(colors) <- tsneMeans$name
    
    if(layout == "tsne" & loop == TRUE) {
        plot <- .plotTsneLoop(graphDF, tsneMeans, colors, d)
    }
    
    if(layout == "tsne" & loop == FALSE) {
        plot <- .plotTsne(graphDF, tsneMeans, colors, d)
    }
    
    if(layout == "igraph" & loop == TRUE) {
        plot <- .plotIgraphLoop(graphDF, colors)
    }
    
    if(layout == "igraph" & loop == FALSE) {
        plot <- .plotIgraph(graphDF, colors)
    }
    
    plot
    return(plot)
})

.plotTsne <- function(graphDF, tsneMeans, colors, d) {
    #setup layout
    layout <- createLayout(graphDF, 'manual', node.positions = tsneMeans)
    layout$x <- tsneMeans$x[sapply(layout$name, function(x) which(tsneMeans$name == x))]
    layout$y <- tsneMeans$y[sapply(layout$name, function(x) which(tsneMeans$name == x))]
    
    #plot
    plot <- ggraph(graph = graphDF, layout = 'manual', node.positions = layout)+
        geom_edge_link(edge_colour="black", aes_string(edge_alpha='weight'))+
        scale_colour_manual(values=colors)+
        geom_node_point(data=d, aes_string(colour='name'), alpha=0.25)+
        geom_node_point(aes_string(colour='name'), size=5)
    
    return(plot)
}

.plotTsneLoop <- function(graphDF, tsneMeans, colors, d) {
    #setup layout
    layout <- createLayout(graphDF, 'manual', node.positions = tsneMeans)
    layout$x <- tsneMeans$x[sapply(layout$name, function(x) which(tsneMeans$name == x))]
    layout$y <- tsneMeans$y[sapply(layout$name, function(x) which(tsneMeans$name == x))]
    
    #plot
    plot <- ggraph(graph = graphDF, layout = 'manual', node.positions = layout)+
        geom_edge_link(edge_colour="black", aes_string(edge_alpha='weight'))+
        geom_edge_loop(
            edge_colour="black",
            aes_string(
                angle=90,
                direction=270,
                strength=50,
                edge_alpha='weight'
            )
        )+
        scale_colour_manual(values=colors)+
        geom_node_point(data=d, aes_string(colour='name'), alpha=0.25)+
        geom_node_point(aes_string(colour='name'), size=5)
    
    return(plot)
}

.plotIgraph <- function(graphDF, colors) {
    plot <- ggraph(graph = graphDF, layout = 'igraph', algorithm = 'kk')+
        geom_edge_link(edge_colour="black", aes_string(edge_alpha='weight'))+
        geom_node_point(aes_string(colour='name'), size=4)+
        scale_colour_manual(values=colors)
    return(plot)
}

.plotIgraphLoop <- function(graphDF, colors) {
    plot <- ggraph(graph = graphDF, layout = 'igraph', algorithm = 'kk')+
        geom_edge_link(edge_colour="black", aes_string(edge_alpha='weight'))+
        geom_node_point(aes_string(colour='name'), size=4)+
        scale_colour_manual(values=colors)+
        geom_edge_loop(
            edge_colour="black",
            aes_string(
                angle=90,
                direction=270,
                strength=1,
                edge_alpha='weight'
            )
        )
}

.swarmTsneMeansProcess <- function(x) {
    spUnsupervised <- getData(x, "spUnsupervised")
    tsneMeans <- getData(spUnsupervised, "tsneMeans")
    colnames(tsneMeans) <- c("name", "x", "y")
    tsneMeans$name <- as.character(tsneMeans$name)
    return(tsneMeans)
}

.swarmTsneProcess <- function(x) {
    spUnsupervised <- getData(x, "spUnsupervised")
    tsne <- getData(spUnsupervised, "tsne")
    classification <- getData(spUnsupervised, "classification")
    
    d <- cbind(as.data.frame(
        tsne[ ,1:2]),
        classification=classification,
        stringsAsFactors=FALSE
    )
    
    colnames(d) <- c("x", "y", "name")
    return(d)
}

.swarmNetworkDF <- function(x) {
    names <- colnames(x)[c(-1,-2)]
    for(o in 1:nrow(x)) {
        ind <- which(x[o, c(-1,-2)] != 0)
        
        if(length(ind) == 1) {
            combs <-  data.frame(V1=names[ind], V2=names[ind], stringsAsFactors=FALSE)
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

.swarmCalculateWeights <- function(graph) {
    t <- as.data.frame(table(graph))
    i <- sapply(t, is.factor)
    t[i] <- lapply(t[i], as.character)
    
    combs <- combn(unique(as.character(t$from)), 2)
    self <- unique(c(t[t$from == t$to, "from"], t[t$from == t$to, "to"]))
    self <- matrix(rep(self, 2), nrow=2, byrow=TRUE)
    combs <- cbind(combs, self)
    
    weight <- c()
    for(i in 1:ncol(combs)) {
        currCombo <- combs[,i]
        if(currCombo[1] == currCombo[2]) {
            w <- as.numeric(subset(t, from %in% currCombo & to %in% currCombo, select="Freq"))
        } else {
            tmp <- t[t$from != t$to, ]
            conns <- subset(tmp, from %in% currCombo & to %in% currCombo)
            w <- sum(conns$Freq)
        }
        weight <- c(weight, w)
    }
    
    c <- as.data.frame(t(combs))
    colnames(c) <- c("from", "to")
    c$weight <- weight
    graph <- subset(c, weight != 0)
    return(graph)
}


