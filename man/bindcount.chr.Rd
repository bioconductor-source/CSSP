\name{bindcount.chr}
\alias{bindcount.chr}
\title{Compute the number of reads overlapping the specified positions for a single chromosome.}
\usage{
bindcount.chr(tagdat, bindpos, fragL = 200, whs = 250)
}
\arguments{
  \item{tagdat}{A \link{numeric} vector of the genome
  coordinates for the starting positions of the aligned
  reads, with positive numbers representing the 5' strand
  and negative numbers representing the 3' strand.}

  \item{bindpos}{A \link{numeric} vector of the genome
  coordinates whose numbers of covering tags are computed.}

  \item{fragL}{A \link{numeric} value for the fragment
  length of the sequencing reads. Default: 200.}

  \item{whs}{A \link{numeric} value for the half window
  size around the binding position. All tags overlapping
  this region are counted. Default: 250.}
}
\value{
A \link{numeric} vector of the numbers of reads overlapping
each position corresponding to "bindpos".
}
\description{
Compute the number of reads overlapping the specified
positions for a single chromosome.
}
\examples{
data( tagdat_chip )
data( bindpos )
bindcount.chr( tagdat_chip[[1]], bindpos[[1]], fragL = 100, whs = 300 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

