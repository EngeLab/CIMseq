
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
#' @param x An spCounts, spUnsupervised, or spSwarm object.
#' @param y An spCounts, spUnsupervised, or spSwarm object.
#' @param z An spCounts object.
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
        p <- .countsErccPlot(
            x,
            y
        )
        p
        return(p)
    }
    if( type == "markers" ) {
        p <- .countsMarkersPlot(
            x,
            y,
            markers
        )
        p
        return(p)
    }
})

#get and process data for ercc plot
.countsErccPlotProcess <- function(x, y) {
    
    sampleType <- c(
        rep(
            "Singlet",
            ncol(getData(x, "counts"))
        ),
        rep(
            "Multuplet",
            ncol(getData(y, "counts"))
        )
    )
    
    counts <- cbind(
        getData(x, "counts"),
        getData(y, "counts")
    )
    counts.ercc <- cbind(
        getData(x, "counts.ercc"),
        getData(y, "counts.ercc")
    )
    
    frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))
    
    d <- data.frame(
        sampleType = factor(sampleType, levels=c("Singlet", "Multuplet")),
        frac.ercc=frac.ercc
    )
    
    #calculate cell number
    d$cellNumber <- median(d[d[,1] == "Singlet", "frac.ercc"]) / d$frac.ercc
    return(d)
}

#plot ercc plot
.countsErccPlot <- function(x, y) {
    
    d <- .countsErccPlotProcess(x, y)
    
    convertTocellNumber <- function(x, d) {
        #100*(1-(median(d[d[,1] == "Singlet", "frac.ercc"])/x))
        100*(median(d[d[,1] == "Singlet", "frac.ercc"])/x)
    }
    
    p <- ggplot(d, aes_string(x='sampleType', y='cellNumber'))+
    geom_jitter()+
    labs(
        x="Sample type",
        y="Cell number",
        title="Fraction ERCC in singlets/doublets"
    )+
    theme_few()+
    theme(
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text=element_text(size=13),
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
    )+
    scale_y_continuous(
        sec.axis = sec_axis(
            trans=~convertTocellNumber(., d),
            name="% ERCC",
            breaks=c(0, seq(60,100,10))
        )
    )
    
    return(p)
}

#get and process data for markers plot
.countsMarkersPlotProcess <- function(x, y, markers) {
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
.countsMarkersPlot <- function(x, y, markers) {
    
    d <- .countsMarkersPlotProcess(x, y, markers)
    
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
    y = NULL,
    type,
    markers = NULL,
    plotUncertainty=TRUE,
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
.unsupClusterPlotProcess <- function(x, plotUncertainty) {
    tsne <- getData(x, "tsne")
    classification <- getData(x, "classification")
    
    d <- cbind(as.data.frame(tsne[ ,1:2]), classification=classification)
    
    if(plotUncertainty) {
        d$uncertainty <- getData(x, "uncertainty")
    } else {
        d$uncertainty <- 3
    }
    
    return(d)
}

#plot clusters plot
.unsupClustersPlot <- function(x, plotUncertainty) {
    
    d <- .unsupClusterPlotProcess(x, plotUncertainty)
    colors <- .setColors()
    
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
.unsupMarkerPlotProcess <- function(x, y, markers, plotUncertainty) {
    tsne <- as.data.frame(getData(x, "tsne"))
    counts.log <- getData(y, "counts.log")
    
    #on the next line you will need a check that the markers exist in the data
    markExpress <- t(counts.log[rownames(counts.log) %in% markers, ])
    markExpress <- apply(markExpress, 2, function(x) {
        (x-min(x))/(max(x)-min(x))
    })
    
    names <- rownames(markExpress)
    markExpress <- as.data.frame(markExpress, row.names=1:nrow(markExpress))
    markExpress$sample <- names
    
    if(plotUncertainty) {
        markExpress$uncertainty <- getData(x, "uncertainty")
    } else {
        markExpress$uncertainty <- 0.0001
    }
    df <- cbind(tsne, markExpress)
    m <- melt(df, id.vars=c("V1", "V2", "sample", "uncertainty"))
    
    return(m)
}

#plot markers plot
.unsupMarkersPlot <- function(x, y, markers, plotUncertainty) {
    
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
        p <- p+geom_point(aes_string(size='uncertainty'))
    } else {
       p <- p+geom_point()
    }
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
#' @importFrom igraph graph_from_data_frame delete_edges set_edge_attr
#'    get.edgelist V
#' @importFrom ggforce theme_no_axes
#' @importFrom utils combn
#' @importFrom ggthemes theme_few

setMethod("spPlot", "spSwarm", function(
    x,
    y,
    z=NULL,
    layout="tsne",
    loop=FALSE,
    markers=NULL,
    edge.cutoff=0,
    min.pval=1,
    min.num.edges=0,
    ...
){
    #x should be an spSwarm object
    #y should be an spUnsupervised object
    #z should be an spCounts object with singlets
    
    #get data
    tsneMeans <- getData(y, "tsneMeans")
    d <- .swarmTsneProcess(y)
    d <- d[order(d$classification), ]

    if(!is.null(markers)) {
        d$color <- .col.from.targets(markers, getData(z, "counts.log"))
    }
    
    graph <- spSwarmPoisson(
        x,
        edge.cutoff=edge.cutoff,
        min.pval=1,
        min.num.edges=0
    )
    
    #convert to igraph
    graphDF <- graph_from_data_frame(graph)
    
    #remove edges according to min.num.edges and min.pval
    remove <- which(graph$weight < min.num.edges | graph$pval > min.pval)
    graphDF <- delete_edges(graphDF, remove)
    
    #add nodes not in graph (due to min.pval argument)
    vertexToAdd <- tsneMeans$classification[!tsneMeans[,1] %in% V(graphDF)$name]
    
    if(length(vertexToAdd) > 0) {
        for(i in 1:length(vertexToAdd)) {
            graphDF <- graphDF + vertex(as.character(vertexToAdd[i]))
        }
    }
    
    #name edges
    graphDF <- set_edge_attr(
        graphDF,
        name="name",
        value=paste(
            get.edgelist(graphDF)[,1],
            get.edgelist(graphDF)[,2],
            sep="->"
        )
    )
    
    #setup colors
    colors <- .setColors()[1:length(unique(tsneMeans$classification))]
    names(colors) <- tsneMeans$classification
    
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

.swarmTsneProcess <- function(x) {
    tsne <- getData(x, "tsne")
    classification <- getData(x, "classification")
    
    d <- cbind(as.data.frame(
        tsne[ ,1:2]),
        classification=classification,
        stringsAsFactors=FALSE
    )
    
    colnames(d) <- c("x", "y", "classification")
    return(d)
}

.plotTsne <- function(graphDF, tsneMeans, colors, d) {
    
    #setup layout
    layout <- create_layout(graphDF, 'manual', node.positions = tsneMeans)
    
    #plot
    plot <- ggraph(
        graph = graphDF,
        layout = 'manual',
        node.positions = layout
    )+
    geom_edge_link(
        edge_colour="black",
        aes_string(edge_width ='factor(weight)'),
        edge_alpha=0.3,
        lineend = "round"
    )+
    geom_node_point(
        data=d,
        aes_string(colour='classification'),
        alpha=0.3
    )+
    geom_node_point(
        data=tsneMeans,
        aes_string(colour='classification'),
        size=5
    )+
    scale_colour_manual(
        name="classification",
        values=colors
    )+
    guides(
        colour = guide_legend(title="Classification"),
        edge_width = guide_legend(title="Weight")
    )+
    theme_few()
        
    return(plot)
}

.plotTsneLoop <- function(graphDF, tsneMeans, colors, d) {
    #setup layout
    layout <- create_layout(graphDF, 'manual', node.positions = tsneMeans)
    
    #plot
    plot <- ggraph(
        graph = graphDF,
        layout = 'manual',
        node.positions = layout
    )+
    geom_edge_link(
        edge_colour="black",
        aes_string(edge_width ='factor(weight)'),
        edge_alpha=0.3,
        lineend = "round"
    )+
    geom_edge_loop(
        edge_colour="black",
        edge_alpha=0.3,
        aes_string(
            angle=90,
            direction=270,
            strength=50,
            edge_width='factor(weight)'
        ),
        lineend = "round"
    )+
    geom_node_point(
        data=d,
        aes_string(colour='classification'),
        alpha=0.3
    )+
    geom_node_point(
        data=tsneMeans,
        aes_string(colour='classification'),
        size=5
    )+
    scale_colour_manual(
        name="classification",
        values=colors
    )+
    guides(
        colour = guide_legend(title="Classification"),
        edge_width = guide_legend(title="Weight")
    )+
    theme_few()
    
    return(plot)
}

.plotIgraph <- function(graphDF, colors) {
    plot <- ggraph(
        graph = graphDF,
        layout = 'igraph',
        algorithm = 'kk'
    )+
    geom_edge_link(
        edge_colour="black",
        aes_string(
            edge_width='factor(weight)'
        ),
        edge_alpha=0.3
    )+
    geom_node_point(
        aes_string(colour='name'),
        size=4
    )+
    scale_colour_manual(values=colors)+
    theme_few()+
    guides(
        colour = guide_legend(title="Classification"),
        edge_width = guide_legend(title="Weight")
    )
    
    return(plot)
}

.plotIgraphLoop <- function(graphDF, colors) {
    plot <- ggraph(
        graph = graphDF,
        layout = 'igraph',
        algorithm = 'kk'
    )+
    geom_edge_link(
        edge_colour="black",
        edge_alpha=0.3,
        aes_string(edge_width='factor(weight)')
    )+
    geom_edge_loop(
        edge_colour="black",
        edge_alpha=0.3,
        aes_string(
            angle=90,
            direction=270,
            strength=1,
            edge_width='factor(weight)'
        )
    )+
    geom_node_point(
        aes_string(colour='name'),
        size=4
    )+
    scale_colour_manual(values=colors)+
    theme_few()+
    guides(
        colour = guide_legend(title="Classification"),
        edge_width = guide_legend(title="Weight")
    )+
    labs(
        x=
    )
}

.col.from.targets <- function(targets, values) {
    
    values <- values[rownames(values) %in% targets, ]
    pal <- .setColors()
    targets <- pal[1:length(targets)]
    targets <- targets[1:(dim(values)[1])]
    if(is.matrix(values)) {
        
        #normalize each genes values so that all genes aer on the same scale
        v <- t(apply(values, 1, function(x) {(x-min(x))/(max(x)-min(x))}))
        
        #calculates the fraction that each gene is expressed in a sample relative to all input genes expression
        fractions <- apply(v, 2, function(x) {x/sum(x)})
        
        #this fixes the Nan problem, due to division with 0, but results in positive expression values for genes samples that had no expression... Do these end up being white in the end?
        fractions[is.nan(fractions)] <- 1.0/(dim(fractions)[1])
        
        #changes the color values from hex to rgb
        targets.rgb <- col2rgb(targets)
        
        #runs a for loop with one iteration per sample
        res <- vector("character", length=length(targets))
        for(i in 1:ncol(values)) {
            
            #for each gene in the current sample, multiply the red, green, and blue values (corresponding to the current gene) with the normalized expression value for that gene. Subtract this from 255 to keep the color within the 0-255 range.
            mytarget.rgb <- 255-t(
                apply(
                    targets.rgb,
                    1,
                    function(x) {
                        (255-x) * v[,i]
                    }
                )
            )
            
            #multiply the rgb color values for this gene by the fraction by which this gene is expressed compared to other target genes
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
            mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(values) {(255-values) * v[i]}))
            #            mytarget.rgb <- rowSums(t(apply(mytarget.rgb, 1, function(values) {values * fractions[,i]})))
            res[i] <- rgb(red=mytarget.rgb[1]/256, green=mytarget.rgb[2]/256, blue=mytarget.rgb[3]/256)
        }
        return(res)
    }
}



