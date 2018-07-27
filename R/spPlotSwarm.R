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
#' @param spSwarm spSwarm; An spSwarm object.
#' @param spUnsupervised spUnsupervised; An spUnsupervised object.
#' @param ... additional arguments to pass on.
#' @return The spPlot function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords plotSwarmGraph
#' @examples
#'
#' #use demo data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#' uObj <- testUns
#' sObj <- testSwa
#'
#' #plot
#' p <- plotSwarmGraph(sObj, uObj)
NULL

#' @rdname plotSwarmGraph
#' @export

setGeneric("plotSwarmGraph", function(
  spSwarm,
  spUnsupervised,
  ...
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

setMethod("plotSwarmGraph", c("spSwarm", "spUnsupervised"), function(
  spSwarm,
  spUnsupervised,
  ...
){
  
  #move data to graph
  spSwarmPoisson(spSwarm, edge.cutoff = 0) %>%
  unite('connection', .data$from, .data$to, sep = "-", remove = FALSE) %>%
  select(.data$from, .data$to, .data$connection, .data$weight, .data$pval) %>%
  graph_from_data_frame(directed = FALSE) %>%
  as_tbl_graph() %>% #better if this could be done directly; avoids importing igraph
  activate(nodes) %>%
  rename('Class' = .data$name) %>%
  #remove edges with 0 weight and coerce to factor
  activate(edges) %>%
  filter(.data$weight > 0) %>%
  mutate('weight' = parse_factor(
    .data$weight,
    levels = unique(.data$weight)
  )) %>%
  #plot base plot with custom layout according to tsne
  ggraph(
    layout = 'manual',
    node.positions = create_layout(
      ., 'manual',
      node.positions = getData(spUnsupervised, "tsneMeans")
    )
  ) +
  #plot edges
  geom_edge_link(
    edge_colour = "black",
    aes_string(edge_width = 'weight'),
    edge_alpha = 0.3,
    lineend = "round"
  ) +
  # add all cells
  geom_node_point(
    data = tidyUnsupervised(spUnsupervised),
    aes_string(
      x = '`t-SNE dim 1`', y = '`t-SNE dim 2`', colour = 'Classification'
    ),
    alpha = 0.3
  ) +
  #add mean point for each class
  geom_node_point(
    data = getData(spUnsupervised, "tsneMeans"),
    aes_string(colour = 'classification'),
    size = 5
  ) +
  #change the color scale
  scale_colour_manual(
    name = "classification",
    values = col40()
  ) +
  theme_few() +
  theme(legend.position = "top", legend.title.align = 0.5) +
  guides(
    colour = guide_legend(title = "Classification", title.position = "top"),
    edge_width = guide_legend(title = "Weight", title.position = "top")
  )
})

#' plotSwarmBarBase
#'
#' Not exported. Provides the base plot for all spSwarm bar type plots.
#'
#' @name plotSwarmBarBase
#' @rdname plotSwarmBarBase
#' @aliases plotSwarmBarBase
#' @param spSwarm spSwarm; An spSwarm object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object with the base data.
#' @author Jason T. Serviss
#' @keywords plotSwarmBarBase
NULL

#' @rdname plotSwarmBarBase

setGeneric("plotSwarmBarBase", function(
    spSwarm,
    ...
){
    standardGeneric("plotSwarmBarBase")
})

#' @rdname plotSwarmBarBase
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom dplyr "%>%" inner_join bind_rows distinct
#' @importFrom ggthemes theme_few

setMethod("plotSwarmBarBase", "spSwarm", function(
    spSwarm,
    ...
){
  #since all of the bar plots show all cell types vs all other cell types and
  #spSwarmPoisson dosen't include duplicate connections, we refomat the data
  #before plotting to include these.
  
  types <- colnames(getData(spSwarm, "spSwarm"))
  from <- rep(types, each = length(types))
  to <- rep(types, length(types))
  d <- tibble(from = from, to = to)
  
  results <- spSwarmPoisson(spSwarm = spSwarm, edge.cutoff = 0)
  
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
#' @aliases plotSwarmEdgeBar
#' @param spSwarm spSwarm; An spSwarm object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmEdgeBar
NULL

#' @rdname plotSwarmEdgeBar

setGeneric("plotSwarmEdgeBar", function(
    spSwarm,
    ...
){
    standardGeneric("plotSwarmEdgeBar")
})

#' @rdname plotSwarmEdgeBar
#' @export
#' @import ggplot2

setMethod("plotSwarmEdgeBar", "spSwarm", function(
    spSwarm,
    ...
){
  plotSwarmBarBase(spSwarm) +
  geom_bar(
    aes_string(x = 'to', y = 'weight', fill = 'to'),
    stat = "identity",
    position = position_dodge(width = 1),
    show.legend = FALSE
  ) +
  geom_label(
  aes_string(x = 'to', y = 'weight + 5', label = 'round(pval, digits = 2)'), #the y axis should be dynamically adjusted
    label.size = 0,
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
#' @param spSwarm spSwarm; An spSwarm object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmPbar
NULL

#' @rdname plotSwarmPbar

setGeneric("plotSwarmPbar", function(
    spSwarm,
    ...
){
    standardGeneric("plotSwarmPbar")
})

#' @rdname plotSwarmPbar
#' @export
#' @import ggplot2

setMethod("plotSwarmPbar", "spSwarm", function(
    spSwarm,
    ...
){
  plotSwarmBarBase(spSwarm) +
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
#' @param spSwarm spSwarm; An spSwarm object.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmHeat
NULL

#' @rdname plotSwarmHeat

setGeneric("plotSwarmHeat", function(
    spSwarm,
    ...
){
    standardGeneric("plotSwarmHeat")
})

#' @rdname plotSwarmHeat
#' @export
#' @import ggplot2
#' @importFrom dplyr desc

setMethod("plotSwarmHeat", "spSwarm", function(
    spSwarm,
    ...
){
  plotSwarmBarBase(spSwarm) +
  geom_tile(aes_string(x = 'from', y = 'to', fill = 'weight')) +
  geom_text(
    aes_string(x = 'from', y = 'to', label = 'round(pval, digits = 2)'),
    colour = "white", size = 5
  ) +
  scale_fill_viridis(option = "E") +
  theme(
    legend.position = "top",
    legend.title.align = 0.5,
    axis.text  = element_blank(),
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
#' @param spSwarm spSwarm; An spSwarm object.
#' @param spCountsMul; An spCounts object containing multiplets.
#' @param genes Character; Genes to be plotted. Can not exceed 20.
#' @param multiplets Character; Multiplets to be plotted.
#' @param fractions Can either be numeric vector, in which case the fractions
#'  are applied to all of the samples, a data.frame mimicking the structure of
#'  the spSwarm slot in an spSwarm object, or NULL, in which case the fractions
#'  are extracted from the spSwarm argument.
#' @param freq Numeric, Length 1 vector indicating the frequency the cost should
#'  be calculated along x.
#' @param ... additional arguments to pass on.
#' @return A ggplot object.
#' @author Jason T. Serviss
#' @keywords plotSwarmGenes
NULL

#' @rdname plotSwarmGenes

setGeneric("plotSwarmGenes", function(
  spSwarm,
  ...
){
  standardGeneric("plotSwarmGenes")
})

#' @rdname plotSwarmGenes
#' @export
#' @import ggplot2
#' @importFrom dplyr desc

setMethod("plotSwarmGenes", "spSwarm", function(
  spSwarm,
  spCountsMul,
  genes,
  multiplets,
  fractions = NULL,
  freq = 10,
  ...
){
  sm <- getData(spSwarm, "syntheticMultiplets")
  cpm <- getData(spCountsMul, "counts.cpm")
  nSyntheticMultiplets <- getData(spSwarm, "arguments")$nSyntheticMultiplets
  
  if(all(is.na(sm))) {
    mess <- paste0(
      "The syntheticMultiplets slot matrix only contains NA values. ",
      "The 'saveSingletData' arg must be TRUE when running spSwarm for this",
      " plot to work."
    )
    stop(mess)
  }
  if(length(genes) > 20) {
    stop("Plotting more than 20 genes at a time is not possible.")
  }
  if(!all(genes %in% rownames(sm))) {
    idx <- which(!genes %in% rownames(sm))
    mess <- paste0(rownames(sm)[idx], " genes were not found in the data")
    stop(mess)
  }
  if(is.data.frame(fractions)) {
    if(nrow(fractions) != length(multiplets)) {
      mess <- paste0(
        "You provided a data.frame of fractions but the number of rows does not",
        " equal the length of the multiplets arg. "
      )
      warning(mess)
    }
  }
  if(is.data.frame(fractions)) {
    if(!all(multiplets %in% rownames(fractions))) {
      mess <- paste0(
        rownames(fractions)[!rownames(fractions) %in% multiplets],
        " were provided via the multiplets arg but are not in the rownames of",
        " the fractions data.frame"
      )
      stop(mess)
    }
  }
  
  #process and merge the synthetic multiplets and real multiplet data.
  processSynthetic <- function(m, nSyntheticMultiplets) {
    m %>%
    matrixStats::rowSums2() %>%
    setNames(rep(1:nSyntheticMultiplets, length(.) / nSyntheticMultiplets)) %>%
    split(., names(.)) %>%
    map(function(x) {names(x) <- genes; x}) %>%
    reduce(bind_rows)
  }
  
  if(is.data.frame(fractions) | is.null(fractions)) {
    
    if(is.null(fractions)) {fractions <- getData(spSwarm, "spSwarm")}
    
    synthetic <- map(multiplets, function(i) {
      t(t(sm[rownames(sm) %in% genes, ]) * as.numeric(fractions[i, ]))
    }) %>%
    map(processSynthetic, nSyntheticMultiplets) %>%
    map2(., multiplets, function(x, y) {
      x %>% add_column(sample = y)
    }) %>%
    bind_rows() %>%
    gather(gene, syntheticValues, -sample)
    
    data <- cpm[rownames(cpm) %in% genes, multiplets] %>%
    matrix_to_tibble("gene") %>%
    gather(sample, count, -gene) %>%
    inner_join(synthetic, by = c("sample", "gene")) %>%
    nest(syntheticValues, .key = "syntheticValues")
    
  } else if(is.numeric(fractions)) {
  
    synthetic <- processSynthetic(
      t(t(sm[rownames(sm) %in% genes, ]) * fractions),
      nSyntheticMultiplets
    )
    
    data <- cpm[rownames(cpm) %in% genes, multiplets] %>%
    matrix_to_tibble("gene") %>%
    gather(sample, count, -gene) %>%
    mutate(syntheticValues = map(gene, function(g) {
      pull(synthetic, g)
    }))
    
  }
  
  #poisson distribution of each synthetic multiplet value
  poissonDistSM <- data %>%
    unnest() %>%
    mutate(pos = list(seq(-10, max(syntheticValues) + 200, freq))) %>%
    mutate(ind.dpois = map2(syntheticValues, pos, ~dpois(round(.y), .x))) %>%
    unnest() %>%
    group_by(gene, sample) %>%
    mutate(ind.dpois.norm = normalizeVec(ind.dpois)) %>%
    ungroup()
  
  #calculate dpois only for the real synthetic multiplet values to be able to
  #show the real mean cost per gene. Ideally a
  #seperate function would calculate dpois and cost.
  realCost <- data %>%
    unnest() %>%
    mutate(dpois.real = map2_dbl(ceiling(count), ceiling(syntheticValues), ~dpois(.x, .y, FALSE))) %>%
    group_by(gene, sample) %>%
    summarize(mean.dpois = mean(dpois.real)) %>%
    ungroup() %>%
    mutate(mean.log = log10(mean.dpois)) %>%
    mutate(cost.real = if_else(is.infinite(mean.log), 323.0052, mean.log * - 1))

  #calculate the entire cost space (blue line)
  entireCost <- data %>%
    unnest() %>%
    mutate(dpois.x = list(round(seq(0, max(syntheticValues), freq)))) %>%
    unnest() %>%
    mutate(dpois = dpois(dpois.x, syntheticValues)) %>%
    group_by(gene, sample, dpois.x) %>%
    summarize(mean = mean(dpois)) %>%
    mutate(mean.log = log10(mean)) %>%
    mutate(cost = if_else(is.infinite(mean.log), 323.005, mean.log * - 1)) %>%
    ungroup() %>%
    mutate(cost.norm = (cost - min(cost)) / (max(cost) - min(cost)))

  max.cost <- max(entireCost$cost)

  #isolate the real multiplet value
  realMultiplet <- data %>%
    select(gene, count, sample) %>%
    distinct()
  
  data %>%
  unnest() %>%
  #plot
  ggplot() +
  #synthetic multiplet values
  geom_histogram(
    aes(syntheticValues, ..density..),
    binwidth = freq, fill = "black"
  ) +
  facet_grid(gene ~ sample, scales = "free_y") +
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
})

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

