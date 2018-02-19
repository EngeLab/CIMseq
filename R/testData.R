#' Test counts data.
#'
#' @title Simulated counts data used to demonstrate sp.scRNAseq package.
#' @docType data
#' @name testCounts
#' @format matrix
#' \describe{
#'     \item{rownames}{Gene names}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @usage testCounts
#' @return Matrix of counts.
#' @examples
#' data(testData)
#'
NULL

#' Test ercc.
#'
#' @title Simulated ERCC data used to demonstrate sp.scRNAseq package.
#' @docType data
#' @name testErcc
#' @format Matrix
#' \describe{
#'     \item{rownames}{ERCC spike-in names.}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @usage testErcc
#' @return Matrix of ERCC counts.
#' @examples
#' data(testData)
#'
NULL

#' spUnsupervised object.
#'
#' @title spUnsupervised object corresponding to test dataset.
#' @docType data
#' @name testUns
#' @format spUnsupervised
#' \describe{
#'     \item{tsne}{See \code{\link{spUnsupervised}}.}
#'     \item{tsneMeans}{See \code{\link{spUnsupervised}}.}
#'     \item{groupMeans}{See \code{\link{spUnsupervised}}.}
#'     \item{classification}{See \code{\link{spUnsupervised}}.}
#'     \item{uncertainty}{See \code{\link{spUnsupervised}}.}
#'     \item{selectInd}{See \code{\link{spUnsupervised}}.}
#' }
#' @usage testUns
#' @return An spUnsupervised object.
#' @examples
#' data(testData)
#'
NULL

#' spSwarm object.
#'
#' @title spSwarm object corresponding to test dataset.
#' @docType data
#' @name testSwa
#' @format spUnsupervised
#' \describe{
#'     \item{spSwarm}{See \code{\link{spSwarm}}.}
#'     \item{costs}{See \code{\link{spSwarm}}.}
#'     \item{convergence}{See \code{\link{spSwarm}}.}
#'     item{arguments}{See \code{\link{spSwarm}}.}
#' }
#' @usage testSwa
#' @return spSwarm object.
#' @examples
#' data(testData)
#'
NULL

