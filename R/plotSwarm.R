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


.ns_legend <- function(data, nonSigCol) {
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

.obsexp_legend <- function(data, pal) {
  p <- data %>%
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
  p <- data %>%
    ggplot() + 
    geom_point(aes(from, to, colour = frac)) + 
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


#' plotSwarmCircos2
#'
#'
#' @name plotSwarmCircos2
#' @rdname plotSwarmCircos2
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
#' @param pal character; A vector including the colour pallete for the score 
#' colours.
#' @param nonSigCol character; Vector of length 1 indicating the colours for 
#' non-significant connections.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author EngeLab
NULL

#' @rdname plotSwarmCircos2

setGeneric("plotSwarmCircos2", function(
  swarm, ...
){
  standardGeneric("plotSwarmCircos2")
})


#' @rdname plotSwarmCircos2
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

setMethod("plotSwarmCircos2", "CIMseqSwarm", function(
  swarm, singlets, multiplets, classOrder = NULL, connectionClass = NULL, 
  alpha = 0.05, weightCut = 0, expectedWeightCut = 0, label.cex = 1, legend = TRUE, 
  pal = colorRampPalette(c("grey95", viridis::viridis(1)))(120)[30:120],
  nonSigCol = "grey95", h.ratio=0.5, maxCellsPerMultiplet = Inf, depleted=FALSE,
  groups=NULL, multiplet.factor=NA, ...
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
  ps <- calculateEdgeStats(swarm=swarm, singlets=singlets, multiplets=multiplets,
                           maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted, groups=groups,
                           multiplet.factor=multiplet.factor) %>%
      mutate(significant = if_else(
                 pval < alpha & weight > weightCut & expected.edges > expectedWeightCut, TRUE, FALSE
             )) %>%
      mutate(col.idx = score) %>%
      mutate(col.idx = if_else(!significant, NA_real_, score)) %>%
      mutate(idx = as.integer(col.idx/max(col.idx, na.rm=T)*(length(pal)-1)+1)) %>%
      mutate(p.col = pal[idx]) %>%
      mutate(p.col = if_else(!significant, nonSigCol, p.col)) %>%
      arrange(idx)

  
  #calculate edge data and add fraction colours
  edges <- longFormConnections(swarm=swarm, singlets=singlets, multiplets=multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted, multiplet.factor=multiplet.factor) %>%
      group_by(connectionName, class) %>%
      mutate(idx = mean(frac)) %>%
      ungroup() %>%
      arrange(match(from, classOrder), match(to, classOrder)) %>%
      mutate(f.col = viridis(101)[idx/max(idx)*100+1])
  
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
      arrange(.closestCircle(match(from, classOrder), match(to, classOrder), match(class, classOrder), 12), .by_group=TRUE) %>%
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
    l2 <- .obsexp_legend(ps, pal)
    l3 <- .frac_legend(edges)
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
  circos.track(
    ylim = c(0, 1), bg.border = "darkgrey", track.height = 0.05, 
    track.margin = c(0.0001, 0.0001), bg.lwd = 0.3,
    panel.fun = function(x, y) {
      sector.index = CELL_META$sector.index
      m <- filter(data, class == sector.index)
      if(nrow(m) == 0) return(NA)
      for(i in 1:nrow(m)) {
        circos.rect(
          xleft = m$position[i], ybottom = 0,
          xright = m$position[i], ytop = 1,
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
#          border = 1,
          col = unique(pull(conn, p.col)),
          lwd=200/length(unique(data$connectionID)),
          h.ratio=h.ratio
      )
  }
  if(legend) par(op)
  
  circos.clear()
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

.closestCircle2 <- function(from, to, max) {
    pos <- to-from
    pos[pos < 0] <- pos[pos<0]+max
    return(-pos)
}

.ns_legend2<- function(data, nonSigCol) {
  p <- data %>%
    ggplot() +
    geom_bar(aes(from, fill = "n.s.")) +
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

.obsexp_legend2 <- function(data, pal) {
    if(nrow(data) == 0) {
        return(0)
    }
    p <- data %>%
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

.frac_legend2 <- function(data) {
  p <- data %>%
    ggplot() + 
    geom_point(aes(meanFrac, idx, colour = meanFrac)) + 
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
#' @param pal character; A vector including the colour pallete for the score 
#' colours.
#' @param nonSigCol character; Vector of length 1 indicating the colours for 
#' non-significant connections.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author EngeLab
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
#' @importFrom tidyr separate
#' @importFrom graphics par layout
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @import circlize

setMethod("plotSwarmCircos", "CIMseqSwarm", function(
  swarm, singlets, multiplets, classOrder = NULL, 
  alpha = 0.05, weightCut = 0, expectedWeightCut = 0, label.cex = 1, legend = TRUE, 
  pal = colorRampPalette(c("grey95", viridis::viridis(1)))(120)[30:120],
  nonSigCol = "grey95", h.ratio=0.9, maxCellsPerMultiplet = Inf, depleted=FALSE, drawFractions=T, multiplet.factor=NA, ...
  ) {
  pval <- weight <- significant <- score <- idx <- p.col <- from <- to <- NULL
  frac <- connectionID <- super <- connectionName <- position <- nr <- NULL
  colour <- NULL
    
  if(!is.null(classOrder)) {
      # Check that supplied classOrder is conformant
      if(!identical(sort(unique(getData(singlets, "classification"))), sort(classOrder))) stop("error, classOrder and singlet classification do not match!")
  } else {
      classOrder <- sort(unique(getData(singlets, "classification")))
  }
    
    adj <- adjustFractions(singlets=singlets, multiplets=multiplets, swarm=swarm, binary=T, maxCellsPerMultiplet=maxCellsPerMultiplet)
    cooc <- matrix(0,nrow=length(classOrder), ncol=length(classOrder), dimnames=list(classOrder, classOrder))
    for(c1 in classOrder) {
        for(c2 in classOrder) {
            if(c1 != c2) {
                cooc[c1, c2] <- sum(adj[,c1] & adj[,c2])/sum(adj[,c1])
                if(sum(adj[,c1] & adj[,c2]) <= weightCut){
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
    count = colSums(adj)[classOrder],
    height = colSums(adj[rowSums(adj) > 1,classOrder])/max(colSums(adj[rowSums(adj) > 1,]))
  )

    
    cooc <- t(apply(cooc, 1, function(x) {x/sum(x)}))
    cooc[is.na(cooc)] <- 0
  #calculate statitistics and connection colors
    ps <- calculateEdgeStats(swarm=swarm, singlets=singlets, multiplets=multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted, multiplet.factor=multiplet.factor) %>%
        mutate(weightOk = weight > weightCut & expected.edges > expectedWeightCut) %>%
      mutate(significant = if_else(
                 pval < alpha & weight > weightCut & expected.edges > expectedWeightCut, TRUE, FALSE
             )) %>%
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
    size.table[,2][size.table[,2] < 0.1] <- 0.3 # Min size is 0.3
    size.table <- as.matrix(size.table)

    

  #calculate edge data and add fraction colours
    edges <- longFormConnections(swarm, singlets, multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet, depleted=depleted, multiplet.factor=multiplet.factor) %>%
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
      mutate(f.col = viridis(101)[(meanFrac-minFrac)/(maxFrac-minFrac)*(100)+1])
  
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
  filtered.data <- data %>% filter(significant) %>% filter(weight > weightCut)
  if(legend & nrow(filtered.data)) {
                                        #layout
        op <- par(mar = par("mar")/2)
        layout(
            matrix(c(1, 2, 3, 5, 5, 5, 4, 4, 4), nrow = 3, byrow = TRUE), 
            widths = c(1, 1, 1, 1, 1), heights = c(1, 8, 2)
        )
        l1 <- .ns_legend2(data, nonSigCol)
        l2 <- .obsexp_legend2(filtered.data, pal)
        l3 <- .frac_legend2(data %>% filter(weight > weightCut))
        l4 <- .class_legend(colours)
  } else {
      if(legend) {
          message(paste0('No significant edges at p=', alpha))
                                        #layout without obsexp
          op <- par(mar = par("mar")/2)
          layout(
              matrix(c(1, 2, 4, 4, 3, 3), nrow = 3, byrow = TRUE), 
              widths = c(1, 1, 1, 1), heights = c(1, 8, 2)
          )
          l1 <- .ns_legend2(data, nonSigCol)
          l3 <- .frac_legend2(data %>% filter(weight > weightCut))
          l4 <- .class_legend(colours)
      }
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
        x=get.cell.meta.data("xcenter", sector.index = colours$class[i]),
        y = 1.5,
      labels = i,
      sector.index = colours$class[i], 1, col = "black",
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
    
  
  
  circos.clear()
  # Reset plotting layout
  if(legend)  par(op)
  par(mfrow=c(1,1))
})

