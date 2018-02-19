#'@include All-classes.R
NULL

#' namedListToTibble
#'
#' Converts a named list to a long data frame.
#'
#' @name namedListToTibble
#' @rdname namedListToTibble
#' @author Jason T. Serviss
#' @param l List. The list to be converted.
#' @keywords namedListToTibble
#' @examples
#'
#' l <- list(a=LETTERS[1:10], b=letters[1:5])
#' output <- namedListToTibble(l)
#'
#' @export
#' @importFrom tibble tibble

namedListToTibble <- function(l) {
    if (length(names(l)) != length(l)) {
        stop("The list you submitted might not be named.")
    }
    if (!is.null(names(l[[1]]))) {
        ni <- gsub(".*\\.(.*)$", "\\1", names(unlist(l)))
        n <- rep(names(l), lengths(l))
        tibble::tibble(
            names = n,
            inner.names = ni,
            variables = unname(unlist(l))
        )
    } else {
        n <- rep(names(l), lengths(l))
        tibble::tibble(
            names = n,
            variables = unname(unlist(l))
        )
    }
}

#' matrix_to_tibble
#'
#' Converts a matrix to a tibble without removing rownames.
#'
#' @name matrix_to_tibble
#' @rdname matrix_to_tibble
#' @author Jason T. Serviss
#' @param data matrix; The matrix to be converted.
#' @keywords matrix_to_tibble
#' @examples
#'
#' m <- matrix(rnorm(20), ncol = 2, dimnames = list(letters[1:10], LETTERS[1:2]))
#' output <- matrix_to_tibble(m)
#'
#' @export
#' @importFrom tibble as_tibble rownames_to_column

matrix_to_tibble <- function(data, rowname = "rowname") {
  data %>%
  as.data.frame() %>%
  rownames_to_column(var = rowname) %>%
  as_tibble()
}

#' normalizeVec
#'
#' Normalizes a vector (x) using x - min(x) / max(x) - min(x).
#'
#' @name normalizeVec
#' @rdname normalizeVec
#' @author Jason T. Serviss
#' @param vec numeric; A numeric vector.
#' @keywords normalizeVec
#' @examples
#'
#' normalizeVec(rnorm(100))
#'
#' @export

normalizeVec <- function(vec) {
  (vec - min(vec)) / (max(vec) - min(vec))
}

#' tidyUnsupervised
#'
#' Tidy spUnsupervised objects.
#'
#' @name tidyUnsupervised
#' @rdname tidyUnsupervised
#' @author Jason T. Serviss
#' @param spUnsupervised; spUnsupervised An spUnsupervised object.
#' @keywords tidyUnsupervised
#' @examples
#'
#' tidyUnsupervised(testUns)
#'
#' @export
#' @importFrom dplyr mutate rename "%>%"

tidyUnsupervised <- function(spUnsupervised) {
  getData(spUnsupervised, "tsne") %>%
  matrix_to_tibble(.) %>%
  mutate(Classification = getData(spUnsupervised, "classification")) %>%
  mutate(Uncertainty = getData(spUnsupervised, "uncertainty")) %>%
  rename(`t-SNE dim 1` = V1, `t-SNE dim 2` = V2, Sample = rowname)
}

#' tidySwarm
#'
#' Tidy spSwarm objects.
#'
#' @name tidySwarm
#' @rdname tidySwarm
#' @author Jason T. Serviss
#' @param spSwarm; spSwarm An spSwarm object.
#' @keywords tidySwarm
#' @examples
#'
#' tidySwarm(testSwa)
#'
#' @export
#' @importFrom dplyr "%>%" full_join left_join
#' @importFrom tibble tibble

tidySwarm <- function(spSwarm) {
  getData(spSwarm, "spSwarm") %>%
  matrix_to_tibble(., rowname = "Sample") %>%
  mutate(Costs = getData(spSwarm, "costs")) %>%
  mutate(Convergence = getData(spSwarm, "convergence"))
}

#spSwarmPoisson(spSwarm, edge.cutoff = 0) %>%
#  full_join(
#    getMultipletsForEdge(spSwarm, edge.cutoff = 0, edges = .),
#    by = c("from", "to")
#  ) %>%
#  left_join(
#    matrix_to_tibble(getData(spSwarm, "spSwarm"), rowname = "multiplet"),
#    by = "multiplet"
#  ) %>%
#  full_join(
#    tibble(
#      multiplet = rownames(getData(spSwarm, "spSwarm")),
#      cost = getData(spSwarm, "costs")
#    ),
#    by = "multiplet"
#  ) %>%
#  full_join(
#    tibble(
#      multiplet = rownames(getData(spSwarm, "spSwarm")),
#      convergence = getData(spSwarm, "convergence")
#    ),
#    by = "multiplet"
#  ) %>%
#  nest(-(from:pval), .key = "Multiplet data")
