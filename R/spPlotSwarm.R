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

