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

matrix_to_tibble <- function(data) {
  data %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble()
}
