% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psoptim.R
\name{pso.2.0}
\alias{pso.2.0}
\title{pso.2.0}
\usage{
pso.2.0(par, fn, gr = NULL, lower = -1, upper = 1,
  swarmInit = NULL, control = list(), ...)
}
\arguments{
\item{par}{numeric; Vector with length defining the dimensionality of the 
optimization problem. For more information see \code{\link[pso]{psoptim}}.}

\item{fn}{function; A function to be minimized. For more information 
see \code{\link[pso]{psoptim}}.}

\item{gr}{function; A function to return the gradient if local search is 
BFGS. For more information see \code{\link[pso]{psoptim}}.}

\item{lower}{numeric; Lower bounds on the variables.}

\item{upper}{numeric; Lower bounds on the variables.}

\item{swarmInit}{matrix; Initial swarm positions. Optional to enhance speed.}

\item{control}{list; A list of control parameters. For more information 
see \code{\link[pso]{psoptim}}. In addition, final swarm positions can be
returned by specifying return.swarm = TRUE.}

\item{...}{Additional arguments to pass on.}
}
\value{
See \code{\link[pso]{psoptim}}.
}
\description{
Particle swarm optimization acceepting a start state and early stopping 
criteria.
}
\author{
Claus Bendtsen modified by Martin Enge
}
