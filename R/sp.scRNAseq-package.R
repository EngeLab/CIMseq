#' A package with purpose of giving a p-value to separations of clusters.
#'
#' It is often a common wish after clustering procedure to assess how
#' significant a seen separation is. This package investigates clusters of two
#' or more groups by projecting all points to a line or curve. The placement on
#' this line is used for classification of the points and the probability of
#' the seen separation just being because of chance is evaluated using
#' a permutation method.
#'
#' \tabular{ll}{ Package: \tab ClusterSignificance\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2016-02-28\cr License: \tab GPL-3\cr }
#'
#' @name ClusterSignificance-package
#' @aliases ClusterSignificance-package ClusterSignificance
#' @docType package
#' @author Author: Jason T. Serviss, Jesper R. Gadin
#' @references Reference to published application note (work in progress)
#' @keywords package
#'
#' @import methods
#'
#' @importFrom grDevices rgb
#' @importFrom graphics lines
#' @importFrom graphics arrows
#' @importFrom utils combn
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices col2rgb
#' @importFrom grDevices n2mfrow
#' @importFrom graphics par
#' @importFrom graphics abline
#' @importFrom graphics hist
#'
NULL

