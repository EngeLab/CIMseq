#script to "simplify" sp.scRNAseq depenencies. Essentially remoevs all plotting
#capability. 

if(!"desc" %in% rownames(installed.packages())) {
  source("https://install-github.me/r-lib/desc")
}

library(desc)
library(purrr)
library(dplyr)
library(devtools)

desc <- description$new()

#remove depend versions
depends <- "R"
purrr::map(depends, ~desc_set_dep(.x, "Depends"))

#remove plotting imports
imports <- c(
  "methods", "mclust", "Rtsne", "pso", "S4Vectors", "readr", "rlang", 
  "matrixStats", "dplyr", "tibble", "purrr", "tidyr", "Rcpp", "future.apply"
)
currentImps <- dplyr::filter(desc$get_deps(), type == "Imports")
impToRm <- dplyr::filter(currentImps, !package %in% imports)$package
purrr::map(impToRm, ~desc_del_dep(.x))

#remove version requirments
purrr::map(imports, ~desc_set_dep(.x, "Imports"))

#remove suggests
currentSugg <- dplyr::filter(desc$get_deps(), type == "Suggests")$package
purrr::map(currentSugg, ~desc_del_dep(.x))

#remove files including plotting functions
f <- list.files(path = "R", pattern = "Plot", full.names = TRUE)
file.remove(f)

#remove vignette
if(file.exists('vignettes')) unlink('vignettes', TRUE)

#remove tests
if(file.exists('tests')) unlink('tests', TRUE)

#remove revdep
if(file.exists('revdep')) unlink('revdep', TRUE)

#remove revdep
if(file.exists('README.md')) unlink('README.md')
if(file.exists('README.Rmd')) unlink('README.Rmd')

#remake c++ accessory files
if(file.exists('src/calculateCost.o')) file.remove('src/calculateCost.o')
if(file.exists('src/RcppExports.cpp')) file.remove('src/RcppExports.cpp')
if(file.exists('src/RcppExports.o')) file.remove('src/RcppExports.o')
if(file.exists('src/sp.scRNAseq.so')) file.remove('src/sp.scRNAseq.so')

#compile attributes
Rcpp::compileAttributes()

#run devtools::document to update docs
devtools::document()

#install simplified version
devtools::install()
