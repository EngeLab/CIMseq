#' A package for analyzing tissue spatiality in scRNAseq data.
#'
#' Description
#'
#' \tabular{ll}{ Package: \tab sp.scRNAseq\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2016-02-28\cr License: \tab GPL-3\cr }
#'
#' @name sp.scRNAseq-package
#' @aliases sp.scRNAseq-package sp.scRNAseq
#' @docType package
#' @author Author: Jason T. Serviss
#' @references Reference to published application note (work in progress)
#' @keywords package
#'
#' @import methods
#' @import ggplot2
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
