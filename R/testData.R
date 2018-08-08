#' Test counts data.
#'
#' @title Counts data used to demonstrate sp.scRNAseq package.
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
#' data(testCounts)
#'
NULL

#' Test counts meta data.
#'
#' @title Meta data for the testCounts dataset.
#' @docType data
#' @name testMeta
#' @format tibble
#' \describe{
#'     \item{sample}{Sample ID}
#'     \item{cellNumber}{Indicates if samples are singlets or multiplets.}
#'     \item{cellTypes}{Indicates the cell types included in the sample.}
#' }
#' @usage testMeta
#' @return Tibble with meta data.
#' @examples
#' data(testMeta)
#'
NULL

#' spCounts object with singlets.
#'
#' @title spCounts object with singlets corresponding to test dataset.
#' @docType data
#' @name test_spCountsSng
#' @format spCounts
#' \describe{
#'     \item{counts}{See \code{\link{spCounts}}.}
#'     \item{counts.log}{See \code{\link{spCounts}}.}
#'     \item{counts.cpm}{See \code{\link{spCounts}}.}
#'     \item{counts.ercc}{See \code{\link{spCounts}}.}
#' }
#' @usage test_spCountsSng
#' @return An spCounts object.
#' @examples
#' data(test_spCountsSng)
#'
NULL

#' spCounts object with multiplets.
#'
#' @title spCounts object with multiplets corresponding to test dataset.
#' @docType data
#' @name test_spCountsMul
#' @format spCounts
#' \describe{
#'     \item{counts}{See \code{\link{spCounts}}.}
#'     \item{counts.log}{See \code{\link{spCounts}}.}
#'     \item{counts.cpm}{See \code{\link{spCounts}}.}
#'     \item{counts.ercc}{See \code{\link{spCounts}}.}
#' }
#' @usage test_spCountsMul
#' @return An spCounts object.
#' @examples
#' data(test_spCountsMul)
#'
NULL

#' spUnsupervised object.
#'
#' @title spUnsupervised object corresponding to test dataset.
#' @docType data
#' @name test_spUnsupervised
#' @format spUnsupervised
#' \describe{
#'     \item{tsne}{See \code{\link{spUnsupervised}}.}
#'     \item{tsneMeans}{See \code{\link{spUnsupervised}}.}
#'     \item{classification}{See \code{\link{spUnsupervised}}.}
#'     \item{uncertainty}{See \code{\link{spUnsupervised}}.}
#'     \item{selectInd}{See \code{\link{spUnsupervised}}.}
#' }
#' @usage test_spUnsupervised
#' @return An spUnsupervised object.
#' @examples
#' data(test_spUnsupervised)
#'
NULL

#' spSwarm object.
#'
#' @title spSwarm object corresponding to test dataset.
#' @docType data
#' @name test_spSwarm
#' @format spSwarm
#' \describe{
#'  \item{spSwarm}{See \code{\link{spSwarm}}.}
#'  \item{costs}{See \code{\link{spSwarm}}.}
#'  \item{convergence}{See \code{\link{spSwarm}}.}
#'  \item{singletIdx}{See \code{\link{spSwarm}}.}
#'  \item{arguments}{See \code{\link{spSwarm}}.}
#' }
#' @usage test_spSwarm
#' @return spSwarm object.
#' @examples
#' data(test_spSwarm)
#'
NULL

