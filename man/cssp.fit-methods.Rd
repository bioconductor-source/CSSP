\docType{methods}
\name{cssp.fit}
\alias{cssp.fit}
\alias{cssp.fit,BinData-method}
\alias{cssp.fit,data.frame-method}
\alias{cssp.fit,data.frame-method,BinData-method}
\title{Fit the CSSP Model.}
\arguments{
  \item{dat}{A \link{data.frame} or \link{BinData-class}
  object containing bin-level chip, input, M and GC
  information. For the data.frame object, the columns must
  contain "chip", "input", "M". For BinData object, the
  slots must contain "tagCount", "input", "M". If "GC" is
  not provided, model will be fitted without using
  gc-Content scores.}

  \item{method}{A \link{character} indicating the method of
  fitting algorithm to be used. "mde" (Default) - minimum
  distance estimation; "gem" - the generalized EM method.}

  \item{p1}{The \link{numeric} value for the lower bound
  for the p-value region where the p-values are assumed to
  be uniformly distributed. Default: 0.5.}

  \item{p2}{The \link{numeric} value for the upper bound
  for the p-value region where the p-values are assumed to
  be uniformly distributed. Default: 0.99.}

  \item{beta.init}{The \link{numeric} value for the
  initializing the size parameter for the background model
  of the ChIP sample. If "NULL", the size parameter of the
  fitted input sample model is used.}

  \item{e0.init}{The \link{numeric} value for initializing
  parameter e0. Default: 0.9.}

  \item{ngc}{An \link{integer} value for the number of
  knots used in the spline model for the gc covariate.
  Default: 9.}

  \item{nite}{An \link{integer} value for the maximum
  number of iterations taken. Default: 50.}

  \item{tol}{A \link{numeric} value for the tolerance for
  convergence. Default: 1e-3.}

  \item{useGrid}{The \link{logical} value indicating
  whether the gridding method is used. If TRUE, the
  covariate space is grided adaptively. This trims down the
  sample size for fitting the regression model when the
  data contains too many observations, and is suggested for
  genome-wide analysis. Default: FALSE.}

  \item{nsize}{A \link{numeric} value for the number of
  bins to be randomly chosen in estimating the
  normalizatiing parameters. If Null (default), all bins
  are used in normalization. For genome wide analysis,
  nsize=5000 is suggested.}

  \item{ncomp}{A \link{numeric} value for the number of
  signal components.}
}
\value{
  \link{CSSPFit-class} A CSSPFit object.
}
\description{
  Fit the CSSP Model.
}
\examples{
data( bin.data )
cssp.fit( bin.data )
cssp.fit( bin.data, method = "gem" )
data( bindata.chr1 )
cssp.fit( bindata.chr1 )
cssp.fit( bindata.chr1, method = "gem", ngc = 1 )
}
\author{
  Chandler Zuo \email{zuo@stat.wisc.edu}
}

