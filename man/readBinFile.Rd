\name{readBinFile}
\alias{readBinFile}
\title{Read the bin-level text files containing ChIP and input sample counts as well as M and GC scores.}
\usage{
readBinFile(type = c("chip", "input", "M", "GC"), fileName)
}
\arguments{
  \item{type}{A \link{character} vector indicating data
  types to be imported.  This vector can contain "chip"
  (ChIP data), "input" (input data), "M" (mappability
  score), "GC" (GC content score). Default:
  c("chip","input","M","GC").}

  \item{fileName}{A \link{character} vector of file names,
  each of which matches each element of "type".  "type" and
  "fileName". This vector should have the same length with
  "type" and corresponding elements in two vectors should
  appear in the same order.}
}
\value{
A \link{data.frame} of the processed bin files, containing
ChIP, input, M and GC in different columns.
}
\description{
Read the bin-level text files containing ChIP and input
sample counts as well as M and GC scores.
}
\note{
"chip","input" and "M" files are all mandatory. "GC" file
is optional.
}
\examples{
data( bindata.chr1 )
pwd <- getwd()
local({
setwd( tempdir() )
on.exit( setwd( pwd ) )
write.table( bindata.chr1[,c(1,4)], file = "chr1_map.txt", sep = "\\t",
row.names = FALSE, col.names = FALSE )
write.table( bindata.chr1[,c(1,5)], file = "chr1_gc.txt", sep = "\\t",
row.names = FALSE, col.names = FALSE )
write.table( bindata.chr1[,c(1,2)], file = "chr1_chip.txt", sep = "\\t",
row.names = FALSE, col.names = FALSE )
write.table( bindata.chr1[,c(1,3)], file = "chr1_input.txt", sep = "\\t",
row.names = FALSE, col.names = FALSE )
readBinFile( fileName = c("chr1_chip.txt", "chr1_input.txt", "chr1_map.txt",
"chr1_gc.txt" ) )
file.remove( paste( "chr1_", c( "chip", "input", "map", "gc" ), ".txt", sep = "" ) )
})
}
\author{
Chandler Zuo\email{zuo@stat.wisc.edu}
}

