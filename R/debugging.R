#'@include All-classes.R
NULL

#' plotSwarmPosition
#'
#' Subtitle
#'
#' Plots the swarm position in a spSwarm object for each reported iteration
#' (specified by reportRate) when report equals TRUE.
#'
#' @name plotSwarmPosition
#' @rdname plotSwarmPosition
#' @aliases plotSwarmPosition
#' @param swarm A CIMseqSwarm object.
#' @param multiplet The name of the multiplet to plot the swarm positions for.
#' @param interval The iterations to be plotted if all recorded iteration
#'    intervals should not be plotted.
#' @param ... additional arguments to pass on.
#' @return The plotSwarmPosition function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords plotSwarmPosition
#' @examples
#'
#' #no test data where swarm was tracked with report arg.
NULL

#' @rdname plotSwarmPosition
#' @export

setGeneric("plotSwarmPosition", function(
  swarm, ...
){
  standardGeneric("plotSwarmPosition")
})

#' @rdname plotSwarmPosition
#' @export
#' @import ggplot2
#' @importFrom dplyr "%>%" filter mutate bind_rows
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom ggthemes theme_few
#' @importFrom readr parse_factor
#' @importFrom rlang .data

setMethod("plotSwarmPosition", "CIMseqSwarm", function(
    swarm,
    multiplet,
    interval = NULL,
    ...
){
  position <- swarmMemberID <- iteration <- NULL
  
  p <- getData(swarm, "stats") %>%
    filter(sample %in% multiplet) %>%
    unnest() %>%
    gather(class, position, -(sample:swarmMemberID)) %>%
    mutate(iteration = as.character(iteration)) %>%
    mutate(
      iteration = parse_factor(
        .data$iteration,
        levels = unique(.data$iteration)
      )
    ) %>%
    mutate(class = factor(.data$class)) %>%
    mutate(position = log2(.data$position + 1)) %>%
    #filter(.data$iteration %in% interval) %>%
    ggplot(aes(position)) +
      geom_density(aes(fill = class, colour = class), alpha = 0.5) +
      facet_grid(sample ~ iteration, scales = "free") +
      scale_fill_manual(values = col40()) +
      scale_colour_manual(values = col40()) +
      ggthemes::theme_few() +
      theme(
        legend.position = "top",
        strip.text.y = element_text(angle = 0)
      ) +
      guides(
        fill = guide_legend(
          ncol = 8, title = "Class", title.position = "top", title.hjust = 0.5
        ),
        colour = FALSE
      ) +
      labs(x = "Swarm position", y = "Density")
  p
  return(p)
})
