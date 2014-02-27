\docType{methods}
\name{qBBT}
\alias{qBBT}
\alias{qBBT,CSSPFit-method}
\title{Compute the quantile estimate for the bin-level poisson parameters.}
\usage{
qBBT(fit, prob, depth = fit@lambday, lower = FALSE)

\S4method{qBBT}{CSSPFit}(fit, prob, depth = fit@lambday, lower = FALSE)
}
\arguments{
  \item{fit}{A \link{CSSPFit-class} object for the CSSP
  model.}

  \item{prob}{A \link{numeric} value for the percentile
  level of bin-level poisson parameters.}

  \item{depth}{A \link{numeric} value for the sequencing
  depth at which the quantile is evaluated.}

  \item{lower}{A \link{logical} value. If TRUE, the lower
  quantile is computed. If FALSE (Default), the upper
  quantile is computed.}
}
\value{
A \link{numeric} value for the percentile of bin-level
poisson parameters.
}
\description{
Compute the quantile estimate for the bin-level poisson
parameters.
}
\examples{
data( sampleFit )
qBBT( sampleFit, prob = 0.99, depth = sampleFit@lambday*0.1 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

