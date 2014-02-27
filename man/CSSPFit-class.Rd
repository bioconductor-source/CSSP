\docType{class}
\name{CSSPFit-class}
\alias{CSSPFit-class}
\title{An S-4 class containing the model fit information for CSSP model.}
\description{
\describe{ \item{lambdax}{Sequencing depth of the input
sample.} \item{lambday}{Sequencing depth of the ChIP
sample.} \item{e0}{The normalization parameter for the ChIP
sample.} \item{pi0}{The pi_0 parameter of CSSP model,
denoting the proportion of bins that are enriched.}
\item{mu.chip}{The vector of the estimated hyper means for
the background model of the ChIP sample.}
\item{mu.input}{The vector of the estimated hyper means for
the input sample.} \item{mean.sig}{The vector of the hyper
means for each signal component.} \item{size.sig}{The
vector of the size parameters for each signal component.}
\item{a}{The size parameter of the input sample model.}
\item{b}{The size parameter of the background model for the
ChIP sample.} \item{p.sig}{The vector of the proportions of
enrichment as each signal component across all enrichment
bins.} \item{prob.zero}{The vector of the prior inflated
probability at 0.} \item{post.p.sig}{The matrix for the
posterior probability of each bin being enriched as a
signal component conditioning on the event that the bin is
enriched. Each column corresponds to one signal component.}
\item{post.p.bind}{Posterior probability of each bin being
enriched.} \item{post.p.zero}{Posterior probability of the
inflated probability at 0.} \item{post.shape.sig}{The
matrix for the shape parameters for the posterior gamma
distributions of bin level poisson parameters, conditioning
on the event that the bins are enriched as each signal
component. Each column corresponds to one signal
component.} \item{post.scale.sig}{The matrix for the scale
parameters of the posterior gamma distributions of bin
level poisson parameters, conditioning on the event that
the bins are enriched as each signal component. Each column
corresponds to one signal component.}
\item{post.shape.back}{The shape parameters for the
posterior gamma distributions of bin level poisson
parameters, conditioning on each bin being enriched.}
\item{post.scale.back}{The scale parameters for the
posterior gamma distributions of bin level poisson
parameters, conditioning on each bin being unenriched.}
\item{n}{The number of mappable bins that are fitted by the
model.} \item{k}{The number of signal components.}
\item{map.id}{The indices for the mappable bins that are
fitted by the model.} \item{pvalue}{The continuously
corrected p-values for a subset of ChIP sample bin counts
against the background model.} \item{cum.pval}{The
cumulative distribution for p-values for a subset of ChIP
sample bin counts against the background model.} }
}
\examples{
showClass("CSSPFit")
}

