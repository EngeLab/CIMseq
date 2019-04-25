#' Test counts data.
#'
#' @title Counts data used to demonstrate CIMseq package.
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
NULL
"testCounts"

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
"testMeta"

#' CIMseqSinglets_test a CIMseqSinglets object.
#'
#' @title CIMseqSinglets object corresponding to test dataset.
#' @docType data
#' @name CIMseqSinglets_test
#' @format CIMseqSinglets
#' \describe{
#'     \item{counts}{See \code{\link{CIMseqSinglets}}.}
#'     \item{counts.log}{See \code{\link{CIMseqSinglets}}.}
#'     \item{counts.cpm}{See \code{\link{CIMseqSinglets}}.}
#'     \item{counts.ercc}{See \code{\link{CIMseqSinglets}}.}
#'     \item{dim.red}{See \code{\link{CIMseqSinglets}}.}
#'     \item{classification}{See \code{\link{CIMseqSinglets}}.}
#' }
#' @return An CIMseqSinglets object.
#' @examples
#' data(CIMseqSinglets_test)
NULL
"CIMseqSinglets_test"

#' CIMseqMultiplets_test a CIMseqMultiplets object.
#'
#' @title CIMseqMultiplets object corresponding to test dataset.
#' @docType data
#' @name CIMseqMultiplets_test
#' @format CIMseqMultiplets
#' \describe{
#'     \item{counts}{See \code{\link{CIMseqMultiplets}}.}
#'     \item{counts.log}{See \code{\link{CIMseqMultiplets}}.}
#'     \item{counts.cpm}{See \code{\link{CIMseqMultiplets}}.}
#'     \item{counts.ercc}{See \code{\link{CIMseqMultiplets}}.}
#'     \item{features}{See \code{\link{CIMseqMultiplets}}.}
#' }
#' @return An CIMseqMultiplets object.
#' @examples
#' data(CIMseqMultiplets_test)
#'
NULL
"CIMseqMultiplets_test"

#' CIMseqSwarm_test a CIMseqSwarm object.
#'
#' @title CIMseqSwarm object corresponding to test dataset.
#' @docType data
#' @name CIMseqSwarm_test
#' @format CIMseqSwarm
#' \describe{
#'  \item{fractions}{See \code{\link{CIMseqSwarm}}.}
#'  \item{costs}{See \code{\link{CIMseqSwarm}}.}
#'  \item{convergence}{See \code{\link{CIMseqSwarm}}.}
#'  \item{singletIdx}{See \code{\link{CIMseqSwarm}}.}
#'  \item{arguments}{See \code{\link{CIMseqSwarm}}.}
#' }
#' @return CIMseqSwarm object.
#' @examples
#' data(CIMseqSwarm_test)
#'
NULL
"CIMseqSwarm_test"
