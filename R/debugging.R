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
#' @param x An spSwarm object.
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
    x,
    ...
){
    standardGeneric("plotSwarmPosition")
})

#' @rdname plotSwarmPosition
#' @export
#' @import ggplot2
#' @importFrom dplyr filter mutate bind_rows
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom ggthemes theme_few
#' @importFrom magrittr %>%
#' @importFrom readr parse_factor
#' @importFrom rlang .data

setMethod("plotSwarmPosition", "spSwarm", function(
    x,
    multiplet,
    interval = NULL,
    ...
){
    #get the reported statistics and subset only those multiplets of interest
    stats <- getData(x, "stats")
    classes <- colnames(getData(x, "spSwarm"))
    names(stats) <- rownames(getData(x, "spSwarm"))
    stats <- stats[names(stats) %in% multiplet]
    
    #reformat all the matrixes as tibbles and add the appropriate names
    ll <- stats %>%
        lapply(.data, `[[`, 4) %>%
        lapply(.data %>%
        lapply(.data %>%
        as_tibble() %>%
        setNames(paste("v", 1:ncol(.data), sep = "")) %>%
        mutate(class = .data$classes)
    ))
    
    ll <- lapply(ll, setNames, paste(1:100, sep = ""))
    
    position <- lapply(1:length(ll), function(x) {
        bind_rows(ll[[x]], .id = "iteration")
    })
    
    p <- position %>%
        setNames(names(stats)) %>%
        bind_rows(.data, .id = "multiplet") %>%
        gather(
            "swarmMember",
            "position",
            -.data$iteration,
            -.data$multiplet,
            -.data$class
        ) %>%
        mutate(
            iteration = parse_factor(
                .data$iteration,
                levels = unique(.data$iteration)
            )
        ) %>%
        mutate(class = factor(.data$class)) %>%
        mutate(position = log2(.data$position + 1)) %>%
        filter(.data$iteration %in% interval) %>%
        as_tibble() %>%
        ggplot(aes(position)) +
            geom_density(aes(fill = class, colour = class), alpha = 0.5) +
            facet_grid(multiplet ~ iteration, scale = "free") +
            scale_fill_manual(values = col64()) +
            scale_colour_manual(values = col64()) +
            theme_few() +
            theme(
                legend.position = "top",
                strip.text.y = element_text(angle = 0),
            ) +
            guides(
                fill = guide_legend(
                    ncol = 8,
                    title = "Class",
                    title.position = "top",
                    title.hjust = 0.5
                ),
                colour = FALSE
            ) +
            labs(
                x = "Swarm position",
                y = "Density"
            )
    
    p
    return(p)
})

swarmPositionjoy <- function(multiplets, iterations = NULL) {


}
