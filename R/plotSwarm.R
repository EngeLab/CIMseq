#'@include All-classes.R
NULL

#' plotSwarmGraph
#'
#' Description.
#'
#' Details.
#'
#' @name plotSwarmGraph
#' @rdname plotSwarmGraph
#' @aliases plotSwarmGraph
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords plotSwarmGraph
#' @examples
#'
#' p <- plotSwarmGraph(
#' CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test
#' )
#'
NULL

#' @rdname plotSwarmGraph
#' @export

setGeneric("plotSwarmGraph", function(
  swarm, singlets, multiplets, ...
){
    standardGeneric("plotSwarmGraph")
})

#' @rdname plotSwarmGraph
#' @export
#' @import ggraph
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom dplyr "%>%" select rename
#' @importFrom rlang .data
#' @importFrom tidyr unite

setMethod(
  "plotSwarmGraph", c("CIMseqSwarm", "CIMseqSinglets", "CIMseqMultiplets"), 
  function(
    swarm, singlets, multiplets, ...
){
  tsneMeans <- means.dim.red(
    getData(singlets, "dim.red"), 
    getData(singlets, "classification")
  )
  
  #move data to graph
  p <- calculateEdgeStats(swarm, singlets, multiplets) %>%
    unite('connection', .data$from, .data$to, sep = "-", remove = FALSE) %>%
    select(.data$from, .data$to, .data$connection, .data$weight, .data$pval) %>%
    graph_from_data_frame(directed = FALSE) %>%
    as_tbl_graph() %>% #better if this could be done directly; avoids importing igraph
    activate(nodes) %>%
    rename('Class' = .data$name) %>%
    #remove edges with 0 weight
    activate(edges) %>%
    filter(.data$weight > 0) %>%
    #plot base plot with custom layout according to tsne
    ggraph(
      layout = 'manual',
      node.positions = create_layout(., 'manual', node.positions = tsneMeans)
    ) +
    #plot edges
    geom_edge_link(
      edge_colour = "black", aes_string(edge_width = 'weight'), 
      edge_alpha = 0.3, lineend = "round"
    ) +
    # add all cells
    geom_node_point(
      data = tidySinglets(singlets),
      aes_string(
        x = '`dim.red dim 1`', y = '`dim.red dim 2`', colour = 'classification'
      ),
      alpha = 0.3
    ) +
    #add mean point for each class
    geom_node_point(
      data = tsneMeans, aes_string(colour = 'classification'), size = 5
    ) +
    #change the color scale
    scale_colour_manual(name = "classification", values = col40()) +
    theme_few() +
    theme(legend.position = "top", legend.title.align = 0.5) +
    guides(
      colour = guide_legend(title = "Classification", title.position = "top"),
      edge_width = guide_legend(title = "Weight", title.position = "top")
    )
    
    p
    return(p)
})

#' plotSwarmBarBase
#'
#' Not exported. Provides the base plot for all spSwarm bar type plots.
#'
#' @name plotSwarmBarBase
#' @rdname plotSwarmBarBase
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object with the base data.
#' @author Jason T. Serviss
NULL

#' @rdname plotSwarmBarBase

setGeneric("plotSwarmBarBase", function(
  swarm, singlets, multiplets, ...
){
    standardGeneric("plotSwarmBarBase")
})

#' @rdname plotSwarmBarBase
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom dplyr "%>%" inner_join bind_rows distinct
#' @importFrom ggthemes theme_few

setMethod(
  "plotSwarmBarBase", c("CIMseqSwarm", "CIMseqSinglets", "CIMseqMultiplets"), 
  function(
    swarm, singlets, multiplets, ...
){
  #since all of the bar plots show all cell types vs all other cell types and
  #calculateEdgeStats dosen't include duplicate connections, we refomat the data
  #before plotting to include these.
  #THIS IS ACTUALLY NOT TRUE ANYMORE. ADJUST
  
  types <- colnames(getData(swarm, "fractions"))
  from <- rep(types, each = length(types))
  to <- rep(types, length(types))
  d <- tibble(from = from, to = to)
  
  results <- calculateEdgeStats(swarm, singlets, multiplets)
  
  join1 <- inner_join(d, results, by = c("from", "to"))
  join2 <- inner_join(d, results, by= c("from" = "to", "to" = "from"))
  
  bind_rows(join1, join2) %>%
  distinct() %>%
  ggplot() +
  theme_few()
})

#' plotSwarmEdgeBar
#'
#'
#' @name plotSwarmEdgeBar
#' @rdname plotSwarmEdgeBar
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
NULL

#' @rdname plotSwarmEdgeBar

setGeneric("plotSwarmEdgeBar", function(
  swarm, singlets, multiplets, ...
){
    standardGeneric("plotSwarmEdgeBar")
})

#' @rdname plotSwarmEdgeBar
#' @export
#' @import ggplot2

setMethod(
  "plotSwarmEdgeBar", c("CIMseqSwarm", "CIMseqSinglets", "CIMseqMultiplets"), 
  function(
    swarm, singlets, multiplets, ...
){
  plotSwarmBarBase(swarm, singlets, multiplets) +
  geom_bar(
    aes_string(x = 'to', y = 'weight', fill = 'to'),
    stat = "identity",
    position = position_dodge(width = 1),
    show.legend = FALSE
  ) +
  #the y axis should be dynamically adjusted. Trying 1/10 weight
  geom_label(
    aes_string(x = 'to', y = 'weight + (weight / 10)', label = 'round(pval, digits = 2)'),
    label.size = 0
  ) +
  facet_grid(from ~ to, scales = "free_x") +
  scale_fill_manual(values = col40()) +
  labs(y = "Weight") +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
})

#' plotSwarmPbar
#'
#'
#' @name plotSwarmPbar
#' @rdname plotSwarmPbar
#' @aliases plotSwarmPbar
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmPbar
NULL

#' @rdname plotSwarmPbar

setGeneric("plotSwarmPbar", function(
  swarm, singlets, multiplets, ...
){
    standardGeneric("plotSwarmPbar")
})

#' @rdname plotSwarmPbar
#' @export
#' @import ggplot2

setMethod(
  "plotSwarmPbar", c("CIMseqSwarm", "CIMseqSinglets", "CIMseqMultiplets"), 
  function(
    swarm, singlets, multiplets, ...
){
  plotSwarmBarBase(swarm, singlets, multiplets) +
  geom_bar(
    aes_string(x = 'to', y = '-log10(pval)', fill = 'to'),
    stat = "identity",
    position = position_dodge(width = 1),
    show.legend = FALSE
  ) +
  facet_grid(from ~ to, scales = "free_x") +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey") +
  scale_fill_manual(values = col40()) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
})

#' plotSwarmHeat
#'
#'
#' @name plotSwarmHeat
#' @rdname plotSwarmHeat
#' @aliases plotSwarmHeat
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmHeat
NULL

#' @rdname plotSwarmHeat

setGeneric("plotSwarmHeat", function(
  swarm, singlets, multiplets, ...
){
    standardGeneric("plotSwarmHeat")
})

#' @rdname plotSwarmHeat
#' @export
#' @import ggplot2
#' @importFrom dplyr desc

setMethod(
  "plotSwarmHeat", c("CIMseqSwarm", "CIMseqSinglets", "CIMseqMultiplets"), 
  function(
    swarm, singlets, multiplets, ...
){
  plotSwarmBarBase(swarm, singlets, multiplets) +
  geom_tile(aes_string(x = 'from', y = 'to', fill = 'weight')) +
  geom_text(
    aes_string(x = 'from', y = 'to', label = 'round(pval, digits = 2)'),
    colour = "white", size = 5
  ) +
  scale_fill_viridis(option = "E") +
  theme(
    legend.position = "top",
    legend.title.align = 0.5,
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  guides(fill = guide_colourbar(title = "Weight", title.position = "top"))
})


#' plotSwarmGenes
#'
#'
#' @name plotSwarmGenes
#' @rdname plotSwarmGenes
#' @aliases plotSwarmGenes
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param genes Character; Genes to be plotted. Can not exceed 20.
#' @param multipletsToPlot Character; Multiplets to be plotted.
#' @param freq Numeric, Length 1 vector indicating the frequency the cost should
#'  be calculated along x.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmGenes
NULL

#' @rdname plotSwarmGenes

setGeneric("plotSwarmGenes", function(
  swarm, ...
){
  standardGeneric("plotSwarmGenes")
})

#' @rdname plotSwarmGenes
#' @export
#' @import ggplot2
#' @importFrom dplyr desc
#' @importFrom purrr reduce map2 map2_dbl
#' @importFrom tidyr nest
#' @importFrom stats dpois
#' @importFrom dplyr case_when

setMethod("plotSwarmGenes", "CIMseqSwarm", function(
  swarm, singlets, multiplets, genes, multipletsToPlot, freq = 10, ...
){
  gene <- syntheticMultipletID <- syntheticValues <- pos <- NULL
  ind.dpois.norm <- dpois.x <- cost.norm <- cost.real <- NULL
  
  sm <- appropriateSinglets(
    singlets, getData(swarm, "singletIdx"), 
    1:nrow(getData(singlets, "counts.cpm"))
  )
  rownames(sm) <- str_replace(rownames(sm), "(.*)\\.[0-9]*", "\\1")
  
  cpm <- getData(multiplets, "counts.cpm")
  nSyntheticMultiplets <- getData(swarm, "arguments")$nSyntheticMultiplets
  selectInd <- getData(multiplets, "features")
  fractions <- getData(swarm, "fractions")
  
  if(length(genes) > 10) {
    stop("Plotting more than 10 genes at a time is not possible.")
  }
  if(!all(genes %in% rownames(sm))) {
    idx <- which(!genes %in% rownames(sm))
    mess <- paste0(genes[idx], " not found in the data")
    stop(mess)
  }
  
  #process synthetic data
  synthetic <- map(multipletsToPlot, function(m) {
    f <- as.numeric(fractions[m, ])
    gNames <- unique(rownames(sm)[rownames(sm) %in% rownames(cpm)[selectInd]])
    
    adjustAccordingToFractions(f, sm) %>%
    multipletSums() %>%
    vecToMat(length(selectInd), nSyntheticMultiplets) %>%
    matrix_to_tibble(drop = TRUE) %>%
    add_column(gene = gNames, .before = 1) %>%
    filter(gene %in% genes) %>%
    gather(syntheticMultipletID, syntheticValues, -gene) %>%
    mutate(syntheticMultipletID = str_replace(
      syntheticMultipletID, ".(.*)", "\\1"
    )) %>%
    add_column(sample = m)
  }) %>%
  reduce(bind_rows)
  
  #process real data and bind synthetic
  if(length(multipletsToPlot) == 1) {
    data <- cpm[rownames(cpm) %in% genes, multipletsToPlot] %>%
      as.data.frame() %>%
      setNames(multipletsToPlot) %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      gather(sample, count, -gene) %>%
      inner_join(synthetic, by = c("sample", "gene")) %>%
      nest(syntheticMultipletID, syntheticValues, .key = "syntheticData")
  } else {
    data <- cpm[rownames(cpm) %in% genes, multipletsToPlot] %>%
      matrix_to_tibble("gene") %>%
      gather(sample, count, -gene) %>%
      inner_join(synthetic, by = c("sample", "gene")) %>%
      nest(syntheticMultipletID, syntheticValues, .key = "syntheticData")
  }
  
  
  #poisson distribution of each synthetic multiplet value
  poissonDistSM <- .pd(data, freq)
  realCost <- .rc(data)
  
  #isolate the real multiplet value
  realMultiplet <- data %>%
    select(gene, count, sample) %>%
    distinct()
  
  #calculate the entire cost space (blue line)
  entireCost <- .ec(data, freq)
  max.cost <- max(entireCost$cost)
  
  p <- data %>%
    unnest() %>%
    #plot
    ggplot() +
    #synthetic multiplet values
    geom_rug(aes(syntheticValues)) +
    facet_grid(gene ~ sample, scales = "free") +
    #this just facilitates the histogram legend
    geom_line(
      aes(x = 0, y = 0, linetype = "Synthetic multiplet values  ")
    ) +
    #poisson distribution of individual synthetic multiplets
    geom_line(
      data = poissonDistSM,
      aes(
        pos, ind.dpois.norm,
        group = interaction(gene, sample, syntheticValues),
        linetype = "Distribution  "
      ),
      size = 0.1, colour = "lightgrey"
    ) +
    #costs
    geom_line(
      data = entireCost,
      aes(dpois.x, cost.norm, linetype = "Costs  "),
      colour = "#3366FF"
    ) +
    #real multiplet value
    geom_segment(
      data = realMultiplet,
      aes(
        x = count, xend = count, y = 0, yend = 1.05,
        linetype = "Real multiplet values  "
      ),
      size = 0.5, colour = "red"
    ) +
    #adds mean cost
    geom_text(
      data = realCost,
      aes(
        x = 0, y = 1.1,
        label = paste0("Mean cost: ", round(cost.real, digits = 2)
      )
      ), colour = "gray13", hjust = 0, size = 3.5
    ) +
    theme_few() +
    labs(y = "Density", x = "CPM") +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~. * max.cost, name = "Cost")
    ) +
    scale_linetype_manual(values = c(1, 1, 2, 1)) +
    guides(linetype = guide_legend(
        title = "",
        override.aes = list(
          fill = rep("white", 4),
          colour = c("#3366FF", "lightgrey", "red", "black")
        )
    )) +
    theme(
      legend.position = "top",
      legend.key = element_blank(),
      legend.title = element_blank()
    )
  
  p
  
  return(p)
})

#poisson distribution of each synthetic multiplet value
.pd <- function(data, freq) {
  syntheticValues <- syntheticValues <- pos <- gene <- NULL
  m <- max(map_dbl(data$syntheticData, function(x) max(x$syntheticValues)))
  
  data %>%
    unnest() %>%
    mutate(pos = list(seq(-10, m + 200, freq))) %>%
    unnest() %>%
    mutate(ind.pois = map2(syntheticValues, pos, function(sv, p) {
      dpois(round(p), round(sv))
    })) %>%
    unnest() %>%
    group_by(sample, gene) %>%
    mutate(ind.dpois.norm = case_when(
      all(ind.pois == 0) ~ ind.pois,
      TRUE ~ normalizeVec(ind.pois)
    )) %>%
    ungroup()
}

#check 
#poissonDistSM %>% 
#  ggplot() +
#  geom_line(aes(pos, ind.dpois.norm, group = syntheticMultipletID)) +
#  facet_grid(gene~sample)

#calculate dpois only for the real synthetic multiplet values to be able to
#show the real mean cost per gene.

.rc <- function(data) {
  syntheticData <- NULL
  data %>%
    mutate(cost.real = map2_dbl(count, syntheticData,
      ~costCalc(round(.x), matrix(unlist(.y$syntheticValues), nrow = 1))
    ))
}

#calculate the entire cost space (blue line)
.ec <- function(data, freq) {
  dpois.x <- syntheticValues <- gene <- dpois <- mean.log <- NULL
  m <- max(map_dbl(data$syntheticData, function(x) max(x$syntheticValues)))
  
  data %>%
    unnest() %>%
    mutate(dpois.x = list(round(seq(0, m + 200, freq)))) %>%
    unnest() %>%
    mutate(dpois = dpois(dpois.x, round(syntheticValues))) %>%
    group_by(gene, sample, dpois.x) %>%
    summarize(mean = mean(dpois)) %>%
    mutate(mean.log = log10(mean)) %>%
    mutate(cost = if_else(is.infinite(mean.log), 323.005, mean.log * (-1))) %>%
    ungroup() %>%
    mutate(cost.norm = normalizeVec(cost))
}

#check
#entireCost %>%
#  ggplot() +
#  geom_line(aes(dpois.x, cost.norm), size = 0.5) +
#  facet_grid(gene ~ sample) +
#  geom_segment(data = realMultiplet, aes(
#      x = count, xend = count, y = 0, yend = 1.05
#  ), size = 0.5, colour = "red")

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
    edges = calculateEdgeStats(
      spSwarm, edge.cutoff = edge.cutoff, min.num.edges = min.num.edges,
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

#' plotSwarmCircos
#'
#'
#' @name plotSwarmCircos
#' @rdname plotSwarmCircos
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param classOrder character; Order of the cell types / classes on the circos.
#' @param connectionClass character; Vector of length 1 specifying the class /
#' cell type to plot or NULL to plot all classes / cell types.
#' @param alpha numeric; Vector of length 1 specifying the which p-values are 
#' considered significant.
#' @param weightCut integer; Vector of length 1 specifying weights below which 
#' the p-value should not be calculated.
#' @param label.cex numeric; Vector of length 1 between [0, 1] indicating the
#'  size of the cell type labels.
#' @param legend logical; indicates if the legends should be plotted.
#' @param theoretical.max integer; See \code{\link{estimateCells}}.
#' @param h.ratio numeric; See \code{\link[circlize]{circos.link}}.
#' @param pal character; A vector including the colour pallete for the score 
#' colours.
#' @param nonSigCol character; Vector of length 1 indicating the colours for 
#' non-significant connections.
#' @param classColour character; Colours for the classes. Order should 
#'  correspond to classOrder argument.
#' @param gap.degree numeric; Controls the amount of space between the classes.
#'  See \code{\link[circlize]{circos.par}}.
#' @param clear logical; Should \code{\link[circlize]{circos.clear}} be called
#'   when finished plotting?
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
NULL

#' @rdname plotSwarmCircos

setGeneric("plotSwarmCircos", function(
  swarm, ...
){
  standardGeneric("plotSwarmCircos")
})

#' @rdname plotSwarmCircos
#' @export
#' @import ggplot2
#' @importFrom dplyr mutate if_else group_by ntile ungroup arrange select filter inner_join desc n "%>%" pull
#' @importFrom viridis viridis
#' @importFrom tibble tibble
#' @importFrom tidyr separate replace_na
#' @importFrom graphics par layout
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @import circlize

setMethod("plotSwarmCircos", "CIMseqSwarm", function(
  swarm, singlets, multiplets, classOrder = NULL, connectionClass = NULL, 
  alpha = 0.05, weightCut = 0, label.cex = 1, legend = TRUE, 
  theoretical.max = NULL, h.ratio = 0.5,
  pal = colorRampPalette(c("grey90", viridis::viridis(1)))(120)[20:110],
  nonSigCol = "grey90", classColour = NULL, gap.degree = NULL, clear = TRUE, ...
){
  pval <- weight <- significant <- score <- idx <- p.col <- from <- to <- NULL
  frac <- connectionID <- super <- connectionName <- position <- nr <- NULL
  colour <- NULL
  
  if(!is.null(classOrder)) {
    # Check that supplied classOrder is conformant
    sngClass <- unique(getData(singlets, "classification"))
    if(!identical(sort(sngClass), sort(classOrder))) {
      stop("classOrder and singlet classification do not match!")
    }
  } else {
    classOrder <- unique(getData(singlets, "classification"))
  }
  
  fractions <- getData(swarm, "fractions")
  if(is.null(classOrder)) classOrder <- unique(getData(singlets, "classification"))
  if(is.null(classColour)) classColour <- col40()[1:length(classOrder)]
  colours <- tibble(
    class = classOrder, 
    colour = classColour,
    nr = 1:length(classOrder),
    combined = paste0("(", nr, ") ", class)
  )
  
  #calculate statitistics and connection colors
  ps <- calculateEdgeStats(
    swarm, singlets, multiplets, theoretical.max = theoretical.max
    ) %>%
    mutate(significant = if_else(
      pval < alpha & weight > weightCut, TRUE, FALSE
    )) %>%
    group_by(significant) %>%
    mutate(idx = ntile(score, length(pal))) %>%
    ungroup() %>%
    mutate(p.col = pal[idx]) %>%
    mutate(p.col = if_else(!significant, nonSigCol, p.col)) %>%
    arrange(idx)
  
  #calculate edge data and add fraction colours
  edges <- longFormConnections(
    swarm, singlets, multiplets, theoretical.max = theoretical.max
  ) %>%
    arrange(match(from, classOrder), match(to, classOrder)) %>%
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
    ungroup() %>%
    select(-sub, -super) %>%
    group_by(class) %>%
    arrange(.closestCircle(
      match(from, classOrder), match(to, classOrder), 
      match(class, classOrder), 12
    ), .by_group=TRUE) %>%
    mutate(position = 1:n()) %>%
    ungroup() %>%
    arrange(significant, to, frac) %>%
    as.data.frame()
  
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
    l2 <- .obsexp_legend(data, pal, alpha)
    l3 <- .frac_legend(data)
    l4 <- .class_legend(colours)
  }
  
  #base circos plot
  posdat <- data %>% 
    group_by(class) %>% 
    summarize(start = min(position), end = max(position)) %>% 
    arrange(match(class, colours$class)) %>%
    full_join(tibble(class = classOrder)) %>%
    replace_na(list(start = 1, end = 2)) %>%
    as.data.frame() %>%
    column_to_rownames("class")
  
  if(is.null(gap.degree)) gap.degree <- 200.0 / length(classOrder)
  
  #initialize circos
  circos.par(gap.degree = gap.degree, cell.padding = c(0, 0))
  circos.initialize(factors = as.character(classOrder), xlim = posdat)
  
  #add labels
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data("sector.index")
      xcenter = get.cell.meta.data("xcenter")
      ycenter = get.cell.meta.data("ycenter")
      circos.text(
        x = xcenter, y = ycenter, 
        labels = filter(colours, class == sector.index)$nr, 
        sector.index = sector.index, col = "grey45", facing = "downward", 
        cex = label.cex
      )
  })
  
  #add colour panels
  circos.trackPlotRegion(
    ylim = c(0, 1), bg.col = pull(colours, colour),
    bg.border = NA, track.height = 0.1
  )
  
  #add fractions
  circos.track(
    ylim = c(0, 1), bg.border = "darkgrey", track.height = 0.05, 
    track.margin = c(0.0001, 0.0001), bg.lwd = 0.3,
    panel.fun = function(x, y) {
      sector.index <- CELL_META$sector.index
      m <- filter(data, class == sector.index) %>%
        group_by(to) %>%
        mutate(frac.pos = sort(position)) %>%
        ungroup()
      if(nrow(m) == 0) return(NA)
      for(i in 1:nrow(m)) {
        circos.rect(
          xleft = m$frac.pos[i], ybottom = 0,
          xright = m$frac.pos[i], ytop = 1,
          sector.index = sector.index,
          border = m$f.col[i], col = m$f.col[i]
        )
      }
    })
  
  #add links
  for(i in 1:length(unique(data$connectionID))) {
    conn <- filter(data, connectionID == unique(data$connectionID)[i])
    if(nrow(conn) != 2) stop("error")
    circos.link(
      pull(conn, class)[1], pull(conn, position)[1],
      pull(conn, class)[2], pull(conn, position)[2],
      col = unique(pull(conn, p.col)),
      lwd = 200 / length(unique(data$connectionID)),
      h.ratio = h.ratio
    )
  }
  if(clear) circos.clear()
  if(legend) par(op)
})

.closestCircle <- function(from, to, class, max) {
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

.ns_legend <- function(data, nonSigCol) {
  connectionID <- NULL
  p <- data %>%
    ggplot() +
    geom_bar(aes(connectionID, fill = "n.s.")) +
    guides(fill = guide_legend(
      label.position = "bottom", title.position = "top",
      frame.colour = "black"
    )) +
    scale_fill_manual(values = nonSigCol) +
    theme(
      legend.position = "top",
      legend.margin = margin(t = 0, b = 0, unit = 'cm'),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.box.background = element_rect(colour = "grey20"),
      legend.box.margin = margin(t = 3, r = 3, b = 1, l = 3, unit = "pt")
    )
  l <- g_legend(p)
  draw_legend(l)
}

.obsexp_legend <- function(data, pal, alpha) {
  pval <- score <- NULL
  p <- data %>%
    filter(pval < alpha) %>%
    ggplot() + 
    geom_point(aes(pval, score, colour = score)) + 
    scale_colour_gradientn(colours = c(pal[1], pal[length(pal)])) +
    guides(colour = guide_colorbar(
      title = "Obs. / Exp.", title.position = "top", title.hjust = 0.5
    )) +
    theme(
      legend.position = "top",
      legend.margin = margin(t = 0, b = 0, unit='cm'),
      legend.key = element_blank(),
      legend.box.background = element_blank()
    )
  l <- g_legend(p)
  draw_legend(l)
}

.frac_legend <- function(data) {
  pval <- score <- frac <- NULL
  p <- data %>%
    ggplot() + 
    geom_point(aes(pval, score, colour = frac)) + 
    scale_colour_viridis_c() +
    guides(colour = guide_colorbar(
      title = "Fractions", title.position = "top", title.hjust = 0.5
    )) +
    theme(
      legend.position = "top",
      legend.margin = margin(b = 0, unit='cm')
    )
  l <- g_legend(p)
  draw_legend(l)
}

.class_legend <- function(colours) {
  combined <- colour <- NULL
  p <- colours %>%
    mutate(combined = parse_factor(combined, levels = combined)) %>%
    ggplot() + 
    geom_bar(aes(combined, fill = colour)) + 
    scale_fill_identity(
      guide = "legend", 
      labels = pull(colours, combined), 
      breaks = pull(colours, colour)
    ) +
    guides(fill = guide_legend(title = NULL)) +
    theme(
      legend.position = "bottom", 
      legend.margin = margin(t = 10, r = 0, b = 50, l = 0, unit = "pt")
    )
  l <- g_legend(p)
  draw_legend(l)
}
