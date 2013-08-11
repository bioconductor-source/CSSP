\name{peakcount}
\alias{peakcount}
\title{Compute the number of aligned reads overlapping the specified peak intervals for the whole genome.}
\usage{
  peakcount(chipdat, inputdat, peakpos, fragL = 200)
}
\arguments{
  \item{chipdat}{A \link{list} of the starting positions of
  the ChIP sample aligned reads for each chromosome. The
  sign of each coordinate represents its strand direction,
  with a positive numbers on the 5' strand and a negative
  numbers on the 3' strand.}

  \item{inputdat}{A \link{list} of the starting positions
  of the input sample aligned reads for each chromosome.
  The sign of each coordinate represents its strand
  direction, with a positive numbers on the 5' strand and a
  negative numbers on the 3' strand.}

  \item{peakpos}{A \link{list} containing the genome
  coordinates for each peak interval on each chromosome.
  Each list component is a 2-column matrix containing the
  left and right boundary of the peak intervals on one
  chromosome.}

  \item{fragL}{A \link{numeric} valkue of the fragment
  length of the aligned reads. Default: 200.}
}
\value{
  A \link{list} of the numbers of reads that overlap the
  corresponding peak intervals.
}
\description{
  Compute the number of aligned reads overlapping the
  specified peak intervals for the whole genome.
}
\examples{
data( peakpos )
data( tagdat_input )
data( tagdat_chip )
peakcount( tagdat_chip, tagdat_input, peakpos, fragL = 100 )
}
\author{
  Chandler Zuo \email{zuo@stat.wisc.edu}
}

