\name{tag2bin}
\alias{tag2bin}
\title{Convert the genome coordinates of aligned reads to bin-level counts for all chromosomes.}
\usage{
tag2bin(tagdat, fragL = 200, binS = 200, prob = 1)
}
\arguments{
  \item{tagdat}{A \link{list} of the genome coordinates for
  starting positions of each read, with positive numbers
  representing the 5' strand and negative numbers
  representing the 3' strand. Each list component
  corresponds to a single chromosome.}

  \item{fragL}{A \link{numeric} value for the fragment
  length of reads. Default: 200.}

  \item{binS}{A \link{numeric} value for the bin-size for
  the bin-level counts to be constructed. Default: 200.}

  \item{prob}{A \link{numeric} value for the proportion of
  randomly sampled reads that will be used to create bin
  data. Default: 1 (use all reads).}
}
\value{
A \link{list} of the bin-level counts for each chromosome.
}
\description{
Convert the genome coordinates of aligned reads to
bin-level counts for all chromosomes.
}
\examples{
data( tagdat_chip )
tag2bin( tagdat_chip, fragL = 100, binS = 100 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

