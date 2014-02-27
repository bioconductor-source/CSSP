\docType{methods}
\name{callpeak}
\alias{callpeak}
\alias{callpeak,CSSPFit-method}
\title{Call enriched bins based on the CSSP model.}
\usage{
callpeak(fit, chip, fold = 1.8, min.count = 0, qval = 0.05, method = "",
  depth = fit@lambday)

\S4method{callpeak}{CSSPFit}(fit, chip, fold = 1.8, min.count = 0,
  qval = 0.05, method = "", depth = fit@lambday)
}
\arguments{
  \item{fit}{A \link{CSSPFit-class} object containing the
  fitted CSSP model.}

  \item{chip}{A \link{numeric} vector containing the bin
  counts for the ChIP sample.}

  \item{depth}{A \link{numeric} value for the sequencing
  depth corresponding to the ChIP sample of the "chip"
  argument. If not provided, sequencing depth of "fit" is
  used.}

  \item{fold}{A \link{numeric} value for the fold change
  threshold for peak calling.}

  \item{min.count}{A \link{numeric} value for the minimum
  ChIP count threshold for peak calling.}

  \item{qval}{A \link{numeric} value for the
  false-discovery rate to be controlled. Default: 0.05.}

  \item{method}{A \link{character} value. By default,
  "min.count" is used to threshold the ChIP bin counts. If
  'method=="post"', "min.count" is used to threshold the
  posterior bin-level poisson intensities.}
}
\value{
A \link{numeric} vector of locations for binding bins.
}
\description{
Call enriched bins based on the CSSP model.
}
\examples{
data( sampleFit )
data( bin.data )
callpeak( sampleFit, chip = bin.data@tagCount, fold = 1, min.count = 0 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

