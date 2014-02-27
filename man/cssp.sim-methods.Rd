\docType{methods}
\name{cssp.sim}
\alias{cssp.sim}
\alias{cssp.sim,CSSPFit-method}
\title{Simulate bin binding intensities according to the posterior distributions of the fitted CSSP model.}
\usage{
cssp.sim(fit, x = fit@lambday)

\S4method{cssp.sim}{CSSPFit}(fit, x = fit@lambday)
}
\arguments{
  \item{x}{A \link{numeric} value for the sequencing depth
  of the ChIP sample at which the new binding intensities
  at simulated.}

  \item{fit}{A \link{CSSPFit-class} class object describing
  the CSSP model.}
}
\value{
A \link{list} object containing

\tabular{ll}{ chip \tab A \link{numeric} vector for the
binding intensities for the ChIP sample.\cr

bind \tab A \link{numeric} vector for the simulated binding
regions.\cr

bind.sig \tab A \link{numeric} vector for the signal
component for each bin.\cr }
}
\description{
Simulate bin binding intensities according to the posterior
distributions of the fitted CSSP model.
}
\examples{
data( sampleFit )
cssp.sim( fit = sampleFit, x = sampleFit@lambday*0.1 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

