\name{peakcount.chr}
\alias{peakcount.chr}
\title{Compute the number of aligned reads overlapping peaks for one chromosome.}
\usage{
peakcount.chr(tagdat, peakpos, fragL = 200, unique = FALSE)
}
\arguments{
  \item{tagdat}{A \link{numeric} vector of the genome
  coordinates for the starting positions of aligned reads.
  The signs of coordinates represent their strand
  direction, with positive numbers representing the 5'
  strand and negative numbers representing the 3' strand.}

  \item{peakpos}{A 2-column \link{matrix} matrix containing
  the left and right position of the peaks for one
  chromosome.}

  \item{fragL}{A \link{numeric} value for the fragment
  length of the sequencing reads. Default: 200.}

  \item{unique}{A \link{logical} value for whether only
  reads mapping to unique nucleotide positions are
  counted.}
}
\value{
A \link{numeric} vector of the number of overlapping tags
for all peaks.
}
\description{
Compute the number of aligned reads overlapping peaks for
one chromosome.
}
\examples{
data( peakpos )
data( tagdat_input )
peakcount.chr( tagdat_input[[1]], peakpos[[1]], fragL = 100 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

