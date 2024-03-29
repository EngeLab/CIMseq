plotSwarmCircos2(swarm3down, cObjSng, cObjMul, weightCut=2, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
arbitraryCircos(swarm3down, cObjSng, cObjMul, weightCut=2, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)




arbitraryCircos <- function(
  swarm, singlets, multiplets, connections, classOrder = NULL, connectionClass = NULL, 
  alpha = 0.05, weightCut = 0, expectedWeightCut = 0, label.cex = 1, legend = FALSE, 
  pal = colorRampPalette(c("grey95", viridis::viridis(1)))(120)[30:120],
  nonSigCol = "grey95", h.ratio=0.5, maxCellsPerMultiplet = Inf, depleted=FALSE, ...
){
  pval <- weight <- significant <- score <- idx <- p.col <- from <- to <- NULL
  frac <- connectionID <- super <- connectionName <- position <- nr <- NULL
  colour <- NULL
    
  fractions <- getData(swarm, "fractions")
  if(!is.null(classOrder)) {
      # Check that supplied classOrder is conformant
      if(!identical(sort(unique(getData(singlets, "classification"))), sort(classOrder))) stop("error, classOrder and singlet classification do not match!")
  } else {
      classOrder <- unique(getData(singlets, "classification"))
  }
  colours <- tibble(
    class = classOrder, 
    colour = col40()[1:length(classOrder)],
    nr = 1:length(classOrder),
    combined = paste0("(", nr, ") ", class)
  )
  
  #calculate statitistics and connection colors
  ps <- calculateEdgeStats(swarm, singlets, multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted) %>%
      mutate(significant = if_else(
                 pval < alpha & weight > weightCut & expected.edges > expectedWeightCut, TRUE, FALSE
             )) %>%
      group_by(significant) %>%
      mutate(idx = ntile(score, length(pal))) %>%
      ungroup() %>%
      mutate(p.col = pal[idx]) %>%
      mutate(p.col = if_else(!significant, nonSigCol, p.col)) %>%
      arrange(idx)

  
  #calculate edge data and add fraction colours
  edges <- longFormConnections(swarm, singlets, multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted) %>%
      arrange(match(from, classOrder),  match(to, classOrder)-match(from, classOrder)) %>%
      mutate(idx = rank(frac)) %>%
      mutate(f.col = viridis(max(idx))[idx]) %>%
      select(-idx)

  
  if(is.null(connectionClass)) {
    filtered <- edges
  } else {
    filtered <- filter(edges, from == connectionClass | to == connectionClass)
  }
  
  if(nrow(filtered) == 0) return("none")
  
  #join all data, select directional connections, arrange and add positions
  data <- filtered %>%
      inner_join(ps, by = c("from", "to")) %>%
      separate(connectionID, into = c("super", "sub"), sep = "\\.", remove = FALSE) %>%
      group_by(sub) %>% 
      filter(pval == min(pval)) %>%
      filter(score == max(score)) %>% # Order makes a difference - IMO selecting max enrichment->pval calc is correct.
      ungroup() %>%
      select(-sub, -super) %>% # Remove sub and super cols
      group_by(class) %>%
      arrange(myClosestCircle(match(from, classOrder), match(to, classOrder), match(class, classOrder), 12), .by_group=TRUE) %>%
      mutate(position = 1:n()) %>%
      ungroup() %>%
                                        #    group_by(class, connectionName) %>% # No sense, arrange does not take group into account
      arrange(significant,  desc(pval), frac) %>% # Why is pval needed here?
                                        #    ungroup() %>%
      as.data.frame()

  size.table <- data.frame(rep(1, length(classOrder)), rep(2, length(classOrder)), row.names=classOrder)
  my.tab <- table(edges$class)
  size.table[names(my.tab),2] <- my.tab/2+1
  size.table <- as.matrix(size.table)
#  size.table <- cbind(rep(1, length(table(edges$class))), table(edges$class)/2)
 
  #add legend
  if(legend) {
    #layout
    op <- par(mar = par("mar")/2)
    layout(
      matrix(c(1, 2, 3, 5, 5, 5, 4, 4, 4), nrow = 3, byrow = TRUE), 
      widths = c(1, 1, 1, 1, 1), heights = c(1, 8, 2)
    )

    #create
    l1 <- .ns_legend(data, nonSigCol)
    l2 <- .obsexp_legend(data, pal)
    l3 <- .frac_legend(data)
    l4 <- .class_legend(colours)
  }
  
  #base circos plot
  class.colors <- col40()[1:length(classOrder)]
  names(class.colors) <- classOrder
  
  #if(is.null(classOrder)) classOrder <- unique(c(edges$from, edges$to))
  #circos.par(track.margin = c(0, 0))
  gap.degree <- 200.0/length(classOrder)
  circos.par(gap.degree=gap.degree, cell.padding=c(0,0)) # FIXME: gap.degree might have to be adjusted!
  circos.initialize(factors=as.character(classOrder), xlim=size.table)
#  circos.initialize(factors=classOrder, xlim=size.table)
  # circos.trackPlotRegion(
  #   ylim = c(0, 1), bg.col = class.colors[sort(names(class.colors))],
  #   bg.border = NA, track.height = 0.1
  # )
  circos.trackPlotRegion(
    ylim = c(0, 1), bg.col = pull(colours, colour),
    bg.border = NA, track.height = 0.1
  )
  #add labels
  for(i in 1:nrow(colours)) { #should be ordered by classOrder
    circos.text(
#      x = mean(range(data$position)), y = 0.5, 
      x = 5, y = 1.5, 
      labels = as.character(pull(colours, nr)[i]), 
#      labels = pull(colours, combined)[i], # Ugly, overlaps links, etc.
      sector.index = pull(colours, class)[i], 1, col = "black",
      facing = "downward", cex = label.cex
    )
  }
  
  #add fractions
#  circos.track(
#    ylim = c(0, 1), bg.border = "darkgrey", track.height = 0.05, 
#    track.margin = c(0.0001, 0.0001), bg.lwd = 0.3,
#    panel.fun = function(x, y) {
#      sector.index = CELL_META$sector.index
#      m <- filter(data, class == sector.index)
#      if(nrow(m) == 0) return(NA)
#      for(i in 1:nrow(m)) {
#        circos.rect(
#          xleft = m$position[i], ybottom = 0,
#          xright = m$position[i], ytop = 1,
#          sector.index = sector.index,
#          border = m$f.col[i], col = m$f.col[i]
#        )
#      }
#    })

  for(i in 1:length(connections)) {
      circos.link(connections$from, connections$to, connections$col, lwd=con$lwd, h.ratio=h.ratio)
  }
  
  #add links
#  for(i in 1:length(unique(data$connectionID))) {
#      conn <- filter(data, connectionID == unique(data$connectionID)[i])
#      if(nrow(conn) != 2) stop("error")
#      circos.link(
#          pull(conn, class)[1], pull(conn, position)[1],
#          pull(conn, class)[2], pull(conn, position)[2],
#          col = unique(pull(conn, p.col)),
#          lwd=200/length(unique(data$connectionID)),
#          h.ratio=h.ratio
#      )
#  }
  if(legend) par(op)
  
  circos.clear()
}

myClosestCircle <- function(from, to, class, max) {
    o <- class == to
    to[o] <- from[o]
    from[o] <- class[o]
    if(any(is.na(from))) {
        cat(paste("NA pos: ", from))
    }
    mid <- as.integer(max/2)
    adj <- mid-from
    from <- (from+adj) %% max
    to <- (to+adj) %% max
    pos <- to-from
    pos[pos < 0] <- -mid-pos[pos < 0]
    pos[pos > 0] <- mid-pos[pos > 0]
    return(pos)
}
