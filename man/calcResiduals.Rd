% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIMseqSwarm.R
\name{calcResiduals}
\alias{calcResiduals}
\title{calcResiduals}
\usage{
calcResiduals(singlets, multiplets, swarm, include = NULL,
  fractions = NULL, ...)
}
\arguments{
\item{singlets}{A CIMseqSinglets object.}

\item{multiplets}{A CIMseqMultiplets object.}

\item{swarm}{A CIMseqSwarm object.}

\item{include}{character;  If residuals should only be calculated for a 
subset of the multiplets, include their names here. Default is to calculate 
for all multiplets.}

\item{fractions}{matrix; A matrix of fractions. By default the fractions in 
the CIMseqSwarm object are used.}

\item{...}{additional arguments to pass on}
}
\value{
Residuals (add more description).
}
\description{
Calculates the residuals for each gene and multiplet after deconvolution
based on the CIMseqSwarm results.
}
\author{
Jason T. Serviss
}
