% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIMseqData.R
\docType{methods}
\name{estimateCells}
\alias{estimateCells}
\alias{estimateCells,CIMseqSinglets-method}
\title{estimateCells}
\usage{
\S4method{estimateCells}{CIMseqSinglets}(singlets, multiplets,
  theoretical.max = Inf, warning = TRUE, ...)
}
\arguments{
\item{singlets}{A CIMseqSinglets object.}

\item{multiplets}{A CIMseqMultiplets object.}

\item{theoretical.max}{integer; Cuts the estimated cell number at the 
provided value such that all values > than the theoretical.max become equal
to the theoretical.max.}

\item{warning}{logical; Indicates if a warning should be issued when all ERCC
counts for a sample are equal to 0. If this warning is issued it can be
resolved by setting the ERCC reads for these samples to NA.}

\item{...}{additional arguments to pass on}
}
\value{
A data frame including the fraction of ercc reads and cell counts for
   each sample.
}
\description{
Subtitle
}
\details{
Uses ERCC data to calculate the fraction of ERCC reads in the samples. In
addition, this function utilizes ERCC data to estimate the cell number
in each sample.
}
\examples{

output <- estimateCells(CIMseqSinglets_test, CIMseqMultiplets_test)

}
\author{
Jason T. Serviss
}
