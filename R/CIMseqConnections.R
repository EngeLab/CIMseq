#'@include All-classes.R
NULL

#' CIMseqConnections
#'
#' Subtitle
#'
#' Description
#'
#' @name CIMseqConnections
#' @rdname CIMseqConnections
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param alpha numeric; Vector of length 1 specifying the which p-values are 
#' considered significant.
#' @param weightCut integer; Vector of length 1 specifying weights below which 
#' the p-value should not be calculated (set to non-zero for enrichment).
#' @param expectedWeightCut integer; Vector of length 1 specifying expected weights below which 
#' the p-value should not be calculated (set to non-zero for depletion).
#' @param maxCellsPerMultiplet integer; Vector of length 1 specifying the maximum size of a multiplet
#' @param depleted logical; Set to TRUE to calculate depletion of connections rather than enrichment
#' @param multiplet.factor numeric; Expected average size of multiplets. Used when ERCC controls are not provided.
#' Should be estimated based on microscopy data. If it is not given, and no ERCC controls are provided, it will be estimated from the data (not recommended).
#' @return CIMseqConnections output.
#' @author Martin Enge
#' @examples
#'
#' #use demo data
#'
NULL

#' @rdname CIMseqConnections
#' @export

setGeneric("CIMseqConnections", function(
  singlets, multiplets, swarm, ...
){
  standardGeneric("CIMseqConnections")
})

#' @importFrom dplyr "%>%" bind_rows mutate
#' @rdname CIMseqConnections
#' @export

setMethod("CIMseqConnections", c("CIMseqSinglets", "CIMseqMultiplets", "CIMseqSwarm"), function(
  singlets, multiplets, swarm, alpha = 0.01, weightCut = 10, expectedWeightCut = 0, maxCellsPerMultiplet=4, depleted=F, multiplet.factor=NA
){
    
    adj <- adjustFractions(singlets=singlets, multiplets=multiplets, swarm=swarm, binary=T, maxCellsPerMultiplet=maxCellsPerMultiplet)
    edgeStats <- calculateEdgeStats(swarm=swarm, singlets=singlets, multiplets=multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted, multiplet.factor=multiplet.factor) %>%
        mutate(weightOk = weight > weightCut & expected.edges > expectedWeightCut) %>%
        mutate(significant = if_else(
                   pval < alpha & weight > weightCut & expected.edges > expectedWeightCut, TRUE, FALSE
               ))
    edges <- longFormConnections(swarm, singlets, multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted, multiplet.factor=multiplet.factor)

    return(new(
        "CIMseqConnections",
        conMatr = adj,
        edgeStats = edgeStats,
        edges = edges,
        alpha = alpha,
        weightCut = weightCut,
        expectedWeightCut = expectedWeightCut,
        maxCellsPerMultiplet = maxCellsPerMultiplet,
        depleted = depleted,
        multiplet.factor = as.numeric(multiplet.factor)
    ))
})


#' plotConnectionCircos
#'
#' Subtitle
#'
#' Description
#'
#' @name plotConnectionCircos
#' @rdname plotConnectionCircos
#' @param connections CIMseqConnections; A CIMseqConnections object.
#' @param classOrder character; Order of the cell types / classes on the circos.
#' @param label.cex numeric; Vector of length 1 between [0, 1] indicating the
#'  size of the cell type labels.
#' @param legend logical; indicates if the legends should be plotted.
#' @param pal character; A vector including the colour pallete for the score 
#' colours.
#' @param nonSigCol character; Vector of length 1 indicating the colours for 
#' non-significant connections.
#' @param h.ratio numeric; Vector of length 1 indicating the curvature of links 
#' @param drawFraction logical; Vector of length 1 indicating the colours for 
#' non-significant connections.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Martin Enge
NULL

#' @rdname plotConnectionCircos

setGeneric("plotConnectionCircos", function(
  connections, ...
){
  standardGeneric("plotConnectionCircos")
})


#' @rdname plotConnectionCircos
#' @export
#' @import ggplot2
#' @importFrom dplyr mutate if_else group_by ntile ungroup arrange select filter inner_join desc n "%>%" pull
#' @importFrom viridis viridis
#' @importFrom tibble tibble
#' @importFrom tidyr separate
#' @importFrom graphics par layout
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @import circlize

setMethod("plotConnectionCircos", "CIMseqConnections", function(
  connections, classOrder = NULL, label.cex = 1, legend = TRUE, 
  pal = colorRampPalette(c("grey95", viridis::viridis(1)))(120)[30:120],
  nonSigCol = "grey95", h.ratio=0.5, drawFractions=T, ...
  ) {
  pval <- weight <- significant <- score <- idx <- p.col <- from <- to <- NULL
  frac <- connectionID <- super <- connectionName <- position <- nr <- NULL
  colour <- NULL

  weightCut <- getData(connections, "weightCut")
  conMatr <- getData(connections, "conMatr")
  if(!is.null(classOrder)) {
      # Check that supplied classOrder is conformant
      if(!identical(sort(unique(getData(singlets, "classification"))), sort(classOrder))) stop("error, classOrder and singlet classification do not match!")
  } else {
      classOrder <- colnames(conMatr)
  }
    
    cooc <- matrix(0,nrow=length(classOrder), ncol=length(classOrder), dimnames=list(classOrder, classOrder))
    for(c1 in classOrder) {
        for(c2 in classOrder) {
            if(c1 != c2) {
                cooc[c1, c2] <- sum(conMatr[,c1] & conMatr[,c2])/sum(conMatr[,c1])
                if(sum(conMatr[,c1] & conMatr[,c2]) <= weightCut){
                    cooc[c1, c2] <- 0
                }
            }
        }
    }

  colours <- tibble(
    class = classOrder, 
    colour = col40()[1:length(classOrder)],
    nr = 1:length(classOrder),
    combined = paste0("(", nr, ") ", class),
    count = colSums(conMatr)[classOrder],
    height = colSums(conMatr)[classOrder]/max(colSums(conMatr))
  )

    
    cooc <- t(apply(cooc, 1, function(x) {x/sum(x)}))
    cooc[is.na(cooc)] <- 0
  #calculate statitistics and connection colors
    ps <- getData(connections, "edgeStats") %>%
        group_by(significant) %>%
        mutate(idx = as.integer((score-min(score))/(max(score)-min(score))*(length(pal)-1)+1)) %>%
        ungroup() %>%
        mutate(p.col = pal[idx]) %>%
        mutate(p.col = if_else(!significant, nonSigCol, p.col)) %>%
        arrange(idx)

    # Select the row with highest enrichment (=score). "name" is the undirectional name.
    ps2 <- ps %>%
        group_by(From1=pmin(from, to), To=pmax(from, to)) %>%
        unite("name", c("From1", "To"), sep="--") %>%
        ungroup() %>%
        group_by(name) %>%
        slice(which.min(pval)) %>%
        ungroup()
    
    size.table <- data.frame(rep(0, length(classOrder)), rep(1.0, length(classOrder)), row.names=classOrder)
    for(my.class in classOrder) {
        my.tab <- c(unlist(ps2 %>% filter(from == my.class, weight > weightCut) %>% select(to)), unlist(ps2 %>% filter(to == my.class) %>% select(from)))
        size.table[my.class,2] <- sum(cooc[my.class, my.tab])
    }
    size.table[,2][size.table[,2] < 0.1] <- 0.3
    size.table <- as.matrix(size.table)

    

  #calculate edge data and add fraction colours
    edges <- getData(connections, "edges") %>%
        group_by(connectionName, class) %>%
        mutate(meanFrac = mean(frac), .by_group=T) %>%
        ungroup() %>%
        select(connectionName, class, from, to, meanFrac) %>%
        distinct()  %>%
        arrange(match(from, classOrder), match(to, classOrder))
    
    if(nrow(edges) == 0) return("none")
    
    edges <- edges %>% mutate(OtherClass = if_else(class == from, true=to, false=from))
    
    
    data <- edges %>%
        inner_join(ps2, by = c("class"="from", "OtherClass"="to")) %>%
        group_by(class) %>%
        arrange(.closestCircle2(match(class, classOrder), match(OtherClass, classOrder), length(classOrder)), .by_group=TRUE) %>%
        mutate(position = 1:n()) %>% # Martin testing stuff
        ungroup() %>%
        group_by(significant) %>%
        mutate(maxFrac=max(meanFrac)) %>%
        mutate(minFrac=min(meanFrac)) %>%
        ungroup() %>%
        arrange(significant,  desc(pval)) %>% # Why is pval needed here?
        as.data.frame()

  class.colors <- col40()[1:length(classOrder)]
  names(class.colors) <- classOrder
  

  class.list <- lapply(classOrder, function(x) {
      ps3 <- ps2 %>% filter(from==x | to==x)  %>% mutate(class=x) %>% mutate(OtherClass = if_else(x == from, true=to, false=from)) %>% arrange(.closestCircle2(match(x, classOrder), match(OtherClass, classOrder), length(classOrder)))
      
      ps3$size <- sapply(1:nrow(ps3), function(i) {
          cooc[ps3$class[i], ps3$OtherClass[i]]
      })
      ps3 %>% mutate(PosStart = cumsum(size)-size) %>% mutate(PosEnd = cumsum(size))
  }) # This section is problematic
  names(class.list) <- classOrder
  
  edges.col <- rbind(edges %>% inner_join(ps2, by=c("class"="from", "OtherClass"="to")), edges %>% inner_join(ps2, by=c("class"="to", "OtherClass"="from"))) %>%
      filter(weightOk) %>%
      mutate(maxFrac=max(meanFrac)) %>%
      mutate(minFrac=min(meanFrac)) %>%
      mutate(f.col = pal[-1][(meanFrac-minFrac)/(maxFrac-minFrac)*(length(pal)-2)+1])
  
  class.list <- lapply(class.list, function(x) {
      x$f.col <- sapply(1:nrow(x), function(i) {
          xtmp <- filter(edges.col, class==x$class[i], OtherClass==x$OtherClass[i])$f.col[1]
          if(is.na(xtmp)) {
              xtmp <- nonSigCol
          }
          xtmp
      })
      x
  })
  
  
    #add legend
    if(legend) {
                                        #layout
        op <- par(mar = par("mar")/2)
        layout(
            matrix(c(1, 2, 3, 5, 5, 5, 4, 4, 4), nrow = 3, byrow = TRUE), 
            widths = c(1, 1, 1, 1, 1), heights = c(1, 8, 2)
        )
        
    #create


        l1 <- .ns_legend2(data, nonSigCol)
        l2 <- .obsexp_legend2(data %>% filter(significant) %>% filter(weight > weightCut), pal)
        l3 <- .frac_legend2(data %>% filter(weight > weightCut))
        l4 <- .class_legend(colours)
    }
    
                                        #base circos plot
    gap.degree <- 200.0/length(classOrder)
    circos.par(gap.degree=gap.degree, cell.padding=c(0,0)) # FIXME: gap.degree might have to be adjusted!
    circos.initialize(factors=as.character(classOrder), xlim=size.table)
    circos.track(
        ylim = c(0, 1), bg.col = 'white',
        bg.border = NA, track.height = 0.14
    )
    circos.track(
      ylim = c(0, 1), bg.col = pull(colours, colour),
      track.margin = c(0.001, 0.0001),
      bg.border = NA, track.height = 0.08
    )
    rmin <- get.cell.meta.data("cell.bottom.radius", track.index = 1)
    rmax <- get.cell.meta.data("cell.top.radius", track.index = 1)
    for(cell in classOrder) {
        draw.sector(get.cell.meta.data("cell.start.degree", sector.index = cell),
                    get.cell.meta.data("cell.end.degree", sector.index = cell),
                    rou1 = rmin+(rmax-rmin)*unlist(filter(colours, class==cell) %>% select(height)),
                    rou2 = rmin,
                    col = '#f7efe3ff',
                    border = "grey25",
                    lwd = 1.5
                    )
    }
  #add labels
  for(i in 1:nrow(colours)) { #should be ordered by classOrder
    circos.text(
      x = 5, y = 1.5, 
      labels = as.character(pull(colours, nr)[i]), 
      sector.index = pull(colours, class)[i], 1, col = "black",
      facing = "downward", cex = label.cex
    )
  }

    

    #add fractions
    if(drawFractions) {
        circos.track(
            ylim = c(0, 1), bg.border = "darkgrey", track.height = 0.04, 
            track.margin = c(0.0001, 0.01), bg.lwd = 0.3,
            panel.fun = function(x, y) {
                sector.index = CELL_META$sector.index
                m <- class.list[[sector.index]]
                if(nrow(m) == 0) return(NA)

                for(i in 1:nrow(m)) {
                    xleft <- m$PosStart[i]
                    xright <- m$PosEnd[i]
                    col <- m$f.col[i]
                    circos.rect(
                        xleft=xleft, ybottom=0,
                        xright=xright, ytop=1,
                        sector.index = sector.index,
                        border = col, col = col
                    )
                }
            }
        )
    }
    
    #add links
    ps2 <- ps2 %>% filter(weight > weightCut) %>% arrange(significant, desc(pval))
    for(i in 1:nrow(ps2)) {
        from.p <- ps2$from[i]
        to.p <- ps2$to[i]
        col <- ps2$p.col[i]
        startp <- unlist(class.list[[from.p]] %>% filter(OtherClass == to.p) %>% select(PosStart, PosEnd))
        endp <- unlist(class.list[[to.p]] %>% filter(OtherClass == from.p) %>% select(PosStart, PosEnd))
        circos.link(from.p, startp, to.p, endp, col=col,h.ratio=h.ratio)
    }
    
    if(legend) par(op)
  
    circos.clear()
})

#' @rdname getEdgeStats

setGeneric("getEdgeStats", function(x){
  standardGeneric("getEdgeStats")
})


#' @rdname getEdgeStats
#' @export

setMethod("getEdgeStats", "CIMseqConnections", function(x) {
    cat('es')
    slot(x, "edgeStats")
})

#' @rdname getEdges

setGeneric("getEdges", function(x){
  standardGeneric("getEdges")
})


#' @rdname getEdges
#' @export

setMethod("getEdges", "CIMseqConnections", function(x) {
    slot(x, "edges")
})


#' @rdname coocHeatmap

setGeneric("coocHeatmap", function(x, pal=c("white", "grey"), ...){
    standardGeneric("coocHeatmap")
})


#' @rdname coocHeatmap
#' @export

setMethod("coocHeatmap", "CIMseqConnections", function(x, pal, ...) {
    heatmap(slot(x, "conMatr"), scale='none', col=pal, ...)
})
