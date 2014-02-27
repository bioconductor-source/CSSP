\docType{methods}
\name{cssp.power}
\alias{cssp.power}
\alias{cssp.power,CSSPFit-method}
\title{Compute the weighted average of bin-wise power conditioning on the fold change and minimal ChIP count requirements.}
\usage{
cssp.power(fit, x, ite = 100, fold = 1, min.count = 10, useC = FALSE,
  qval = 0.05)

\S4method{cssp.power}{CSSPFit}(fit, x, ite = 100, fold = 1,
  min.count = 10, useC = FALSE, qval = 0.05)
}
\arguments{
  \item{x}{A \link{numeric} value for the sequencing depth
  of the ChIP sample at which the power is evaluated.}

  \item{fit}{A \link{CSSPFit-class} object for the CSSP
  model.}

  \item{ite}{A \link{integer} value for the number of
  iterations used for Monte-Carlo evaluation.}

  \item{fold}{A \link{numeric} value for the fold change
  threshold.}

  \item{min.count}{A \link{numeric} value for the minimal
  count threshold.}

  \item{useC}{A \link{logical} value. Whether the function
  will be evaluated using C. Default: FALSE.}

  \item{qval}{A \link{numeric} value for the q-value for
  FDR control. Default: 0.05.}
}
\value{
A \link{numeric} value for the weighted average of bin
power conditioning on the minimal count and fold change
thresholds.
}
\description{
Compute the weighted average of bin-wise power conditioning
on the fold change and minimal ChIP count requirements.
}
\examples{
data( sampleFit )
cssp.power( sampleFit, x = sampleFit@lambday*0.1, min.count = 0, fold = 2,
useC = TRUE )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

