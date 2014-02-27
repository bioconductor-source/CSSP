#'@name tag2bin.chr
#'@title Convert the genome coordinates of aligned reads to bin-level counting data for a single chromosome.
#'
#'@param tagdat A \link{numeric} vector of genome coordinates for the starting positions of the aligned reads, with positive numbers representing the 5' strand and negative numbers representing the 3' strand.
#'@param fragL A \link{numeric} value of the fragment length for the reads. Default: 200.
#'@param binS A \link{numeric} value of the bin-size for the bin-level data to be constructed. Default: 200.
#'@return A \link{numeric} vector of the counts for each bin.
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( tagdat_chip )
#'tag2bin.chr( tagdat_chip[[1]], fragL = 100, binS = 100 )
#'@useDynLib CSSP
#'@export
tag2bin.chr <- function(tagdat,fragL=200,binS=200)
{
  if(class(tagdat)=="integer")tagdat <- as.double(tagdat)
  return (.Call("tag2bin_c",tagdat,fragL,binS,max(abs(tagdat))+fragL ,PACKAGE="CSSP"))
}

#' @name tag2bin
#'@title Convert the genome coordinates of aligned reads to bin-level counts for all chromosomes.
#'
#'@param tagdat A \link{list} of the genome coordinates for starting positions of each read, with positive numbers representing the 5' strand and negative numbers representing the 3' strand. Each list component corresponds to a single chromosome.
#'@param fragL A \link{numeric} value for the fragment length of reads. Default: 200.
#'@param binS A \link{numeric} value for the bin-size for the bin-level counts to be constructed. Default: 200.
#'@param prob A \link{numeric} value for the proportion of randomly sampled reads that will be used to create bin data. Default: 1 (use all reads).
#'@return A \link{list} of the bin-level counts for each chromosome.
#'@useDynLib CSSP
#'@examples 
#'data( tagdat_chip )
#'tag2bin( tagdat_chip, fragL = 100, binS = 100 )
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@export
tag2bin <- function(tagdat,fragL=200,binS=200,prob=1)
{
  bindat <- NULL
  fragL <- as.double(fragL)
  binS <- as.double(binS)
  for(i in seq_along(tagdat))
    {
      if(class(tagdat[[i]])=="integer")tagdat[[i]]<-as.double(tagdat[[i]])
      tagdat[[i]]=sample(tagdat[[i]],as.integer(length(tagdat[[i]])*prob))
      bindat[[i]] <- .Call("tag2bin_c",tagdat[[i]],fragL,binS,max(abs(tagdat[[i]]))+fragL ,PACKAGE="CSSP")
    }
  names(bindat) <- names(tagdat)
  return(bindat)
}

#'@name bindcount.chr
#'@title Compute the number of reads overlapping the specified positions for a single chromosome.
#'
#'@param tagdat A \link{numeric} vector of the genome coordinates for the starting positions of the aligned reads, with positive numbers representing the 5' strand and negative numbers representing the 3' strand.
#'@param bindpos A \link{numeric} vector of the genome coordinates whose numbers of covering tags are computed.
#'@param fragL A \link{numeric} value for the fragment length of the sequencing reads. Default: 200.
#'@param whs A \link{numeric} value for the half window size around the binding position. All tags overlapping this region are counted. Default: 250.
#'@return A \link{numeric} vector of the numbers of reads overlapping each position corresponding to "bindpos".
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@useDynLib CSSP
#'@examples
#'data( tagdat_chip )
#'data( bindpos )
#'bindcount.chr( tagdat_chip[[1]], bindpos[[1]], fragL = 100, whs = 300 )
#'@export
bindcount.chr <- function(tagdat,bindpos,fragL=200,whs=250)
{
  return(.Call("bindcount_c",as.double(tagdat),as.double(bindpos),fragL,whs=whs,PACKAGE="CSSP"))
}

#'@name bindcount
#'@title Compute the number of reads overlapping the specified positions for the whole genome.
#'
#'@param chipdat A \link{list} of the starting coordinates for aligned reads for all chromosomes, with positive numbers representing the 5' strand and negative numbers representing the 3' strand.
#'@param inputdat A \link{list} of the starting coordinates for aligned reads for the input sample for all chromosomes, with positive numbers representing the 5' strand and negative numbers representing the 3' strand.
#'@param bindpos A \link{list} of genome coordinates for each chromosome whose numbers of covering tags are computed.
#'@param fragL A \link{numeric} value for the fragment length of the aligned reads. Default: 200.
#'@param whs A \link{numeric} value for the half window size around the binding position. All tags overlapping this region are counted. Default: 250.
#'@return A \link{list} of the number of overlapping tags for all position. Each list is a data.frame corresponding to a single chromosome, containing:
#'
#'\tabular{ll}{
#'chip \tab The number of ChIP sample reads overlapping each position.\cr
#'input \tab The number of input sample reads overlapping each position.\cr
#'}
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( tagdat_input )
#'data( tagdat_chip )
#'data( bindpos )
#'bindcount( tagdat_chip, tagdat_input, bindpos, fragL = 100, whs = 300 )
#'@export
bindcount <- function(chipdat,inputdat,bindpos,fragL=200,whs=250)
  {
    bindCount <- NULL
    for(i in names(bindpos))
      {
        if(!is.null(chipdat[[i]]))
          {
            bindCount[[i]]$chip <- bindcount.chr(chipdat[[i]],bindpos[[i]],fragL=fragL,whs=whs)
          }else{
            bindCount[[i]]$chip <- NULL
          }
        if(!is.null(inputdat[[i]]))
          {
            bindCount[[i]]$input <- bindcount.chr(inputdat[[i]],bindpos[[i]],fragL=fragL,whs=whs)
          }else{
            bindCount[[i]]$input <- NULL
          }
      }
    return (bindCount)
  }


#'@name peakcount
#'@title Compute the number of aligned reads overlapping the specified peak intervals for the whole genome.
#'
#'@param chipdat A \link{list} of the starting positions of the ChIP sample aligned reads for each chromosome. The sign of each coordinate represents its strand direction, with a positive numbers on the 5' strand and a negative numbers on the 3' strand.
#'@param inputdat A \link{list} of the starting positions of the input sample aligned reads for each chromosome. The sign of each coordinate represents its strand direction, with a positive numbers on the 5' strand and a negative numbers on the 3' strand.
#'@param peakpos A \link{list} containing the genome coordinates for each peak interval on each chromosome. Each list component is a 2-column matrix containing the left and right boundary of the peak intervals on one chromosome.
#'@param fragL A \link{numeric} value of the fragment length of the aligned reads. Default: 200.
#'@param unique A \link{logical} value for whether only reads mapping to unique nucleotide positions are counted.
#'@return A \link{list} of the numbers of reads that overlap the corresponding peak intervals. 
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( peakpos )
#'data( tagdat_input )
#'data( tagdat_chip )
#'peakcount( tagdat_chip, tagdat_input, peakpos, fragL = 100 )
#'@export
peakcount <- function(chipdat,inputdat,peakpos,fragL=200,unique=FALSE)
  {
    bindCount <- NULL
    for(i in names(peakpos))
      {
        if(!is.null(chipdat[[i]]))
          {
            bindCount[[i]]$chip <- peakcount.chr(chipdat[[i]],peakpos[[i]],fragL=fragL,unique=unique)
          }else{
            bindCount[[i]]$chip <- NULL
          }
        if(!is.null(inputdat[[i]]))
          {
            bindCount[[i]]$input <- peakcount.chr(inputdat[[i]],peakpos[[i]],fragL=fragL,unique=unique)
          }else{
            bindCount[[i]]$input <- NULL
          }
      }
    return (bindCount)
  }

#'@name peakcount.chr
#'@title Compute the number of aligned reads overlapping peaks for one chromosome.
#'
#'@param tagdat A \link{numeric} vector of the genome coordinates for the starting positions of aligned reads. The signs of coordinates represent their strand direction, with positive numbers representing the 5' strand and negative numbers representing the 3' strand.
#'@param peakpos A 2-column \link{matrix} matrix containing the left and right position of the peaks for one chromosome.
#'@param fragL A \link{numeric} value for the fragment length of the sequencing reads. Default: 200.
#'@param unique A \link{logical} value for whether only reads mapping to unique nucleotide positions are counted.
#'@return A \link{numeric} vector of the number of overlapping tags for all peaks. 
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@useDynLib CSSP
#'@examples 
#'data( peakpos )
#'data( tagdat_input )
#'peakcount.chr( tagdat_input[[1]], peakpos[[1]], fragL = 100 )
#'@export
peakcount.chr <- function(tagdat,peakpos,fragL=200,unique=FALSE)
  {
      if( unique )
          return(.Call("peakcount_uniq",as.numeric(tagdat),as.numeric(peakpos[,1]),as.numeric(peakpos[,2]),fragL,PACKAGE="CSSP"))
      else
          return(.Call("peakcount_c",as.numeric(tagdat),as.numeric(peakpos[,1]),as.numeric(peakpos[,2]),fragL,PACKAGE="CSSP"))
  }

#'@name readBinFile
#'@title Read the bin-level text files containing ChIP and input sample counts as well as M and GC scores.
#'
#'@param type A \link{character} vector indicating data types to be imported.  This vector can contain "chip" (ChIP data), "input" (input data), "M" (mappability score), "GC" (GC content score). Default: c("chip","input","M","GC").
#'@param fileName A \link{character} vector of file names, each of which matches each element of "type".  "type" and "fileName". This vector should have the same length with "type" and corresponding elements in two vectors should appear in the same order.
#'@note "chip","input" and "M" files are all mandatory. "GC" file is optional.
#'@return A \link{data.frame} of the processed bin files, containing ChIP, input, M and GC in different columns.
#'@author Chandler Zuo\email{zuo@@stat.wisc.edu}
#'@examples
#'data( bindata.chr1 )
#'pwd <- getwd()
#'local({
#'setwd( tempdir() )
#'on.exit( setwd( pwd ) )
#'write.table( bindata.chr1[,c(1,4)], file = "chr1_map.txt", sep = "\t",
#' row.names = FALSE, col.names = FALSE )
#'write.table( bindata.chr1[,c(1,5)], file = "chr1_gc.txt", sep = "\t",
#'row.names = FALSE, col.names = FALSE )
#'write.table( bindata.chr1[,c(1,2)], file = "chr1_chip.txt", sep = "\t",
#'row.names = FALSE, col.names = FALSE )
#'write.table( bindata.chr1[,c(1,3)], file = "chr1_input.txt", sep = "\t",
#'row.names = FALSE, col.names = FALSE )
#'readBinFile( fileName = c("chr1_chip.txt", "chr1_input.txt", "chr1_map.txt",
#'"chr1_gc.txt" ) )
#'file.remove( paste( "chr1_", c( "chip", "input", "map", "gc" ), ".txt", sep = "" ) )
#'})
#'@export
readBinFile <- function(type=c("chip","input","M","GC"),fileName)
  {
    existChip <- (length(which(type == "chip")) > 0)
    existInput <- (length(which(type == "input")) > 0)
    existM <- (length(which(type == "M")) > 0)
    existGC <- (length(which(type == "GC")) > 0)
    if(existChip*existInput*existM==0)
      stop("ChIP, input and M files must be provided")
    if(length(fileName)!=length(type))
      stop("fileName and type must have same length")
    input <- read.table(fileName[which(type=="input")])
    input <- input[,ncol(input)-(1:0)]
    chip <- read.table(fileName[which(type=="chip")])
    chip <- chip[,ncol(chip)]
    map <- read.table(fileName[which(type=="M")])
    map <- map[,ncol(map)]
    if(existGC)
      {
        gc <- read.table(fileName[which(type=="GC")])
        gc <- gc[,ncol(gc)]
        n <- min(c(nrow(input),length(chip),length(map),length(gc)))
        outData <- data.frame(cbind(input[1:n,],chip[1:n],map[1:n],gc[1:n]))
        names(outData) <- c("pos","input","chip","M","GC")
      }else{
        n <- min(c(nrow(input),length(chip),length(map)))
        outData <- data.frame(cbind(input[1:n,],chip[1:n],map[1:n]))
        names(outData) <- c("pos","input","chip","M")
      }
    return(outData)
  }


#'@name createBinData
#'@title Create a BinData object by merging lists of ChIP and input bin data with external M and GC text files.
#'
#'@description This function create a BinData object by merging ChIP and input bin-level counts with external M/GC/N text files.
#'@param dat.chip Either a \link{list} of the ChIP bin level data for each chromosome, or a \link{character} string of the file name including the ChIP bin level data. If the ChIP bin level file name is provided, the file must contain at least two columns, where the chromosome information is in the first column, and the bin level counts are in the last column.
#'@param dat.input A \link{list} of the input bin level data for each chromosome, or a \link{character} string for the input bin level data counts. The structure is the same as "dat.chip".
#'@param chrlist A \link{list} of the chromosomes that is imported. If "NULL", all chromosomes specified by "name(dat.chip)" are imported.
#'@param mfile A \link{character} value. If "m.suffix=NULL", this is the file name of the genome-wide M file. Otherwise, this is the common prefix (including relative path) for all chromosome-level M files.
#'@param gcfile A \link{character} value. If "gc.suffix=NULL", this is the file name of the genome-wide GC file. Otherwise, this is the common prefix (including relative path) for all chromosome-level GC files.
#'@param nfile A \link{character} value. If "n.suffix=NULL", this is the file name of the genome-wide N file. Otherwise, this is the common prefix (including relative path) for all chromosome-level N files.
#'@param m.suffix A \link{character} value. If not NULL, this is the suffix of the chromosome-wise M files. The chromosome-level file has to be named "chrX_m.suffix".
#'@param gc.suffix A \link{character} value. If not NULL, this is the suffix of the chromosome-wise GC files. The chromosome-level file has to be named "chrX_gc.suffix".
#'@param n.suffix A \link{character} value. If not NULL, this is the suffix of the chromosome-wise N files. The chromosome-level file has to be named "chrX_n.suffix".
#'@param dataType A \link{character} value of either "unique" or "multi".
#'@note When .suffix is null, the corresponding genome-wise file must have three columns, with the first column being the chromosome names, the second column being the genome coordinates, and the third column being the corresponding scores. In contrast, when .suffix is not null, then each chromosome-level M/GC/N file should only contain two columns, with the first column being the genome coordinates and the second column being the scores.
#'@return A \link{BinData-class} object.
#'@example inst/tests/createBinData-Ex.R
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@export
createBinData <- function(dat.chip,dat.input,mfile,gcfile,nfile,m.suffix=NULL,gc.suffix=NULL,n.suffix=NULL,chrlist=NULL,dataType="unique")  {

  dataFrameToList <- function( indat ) {
    chrlist <- unique( indat[,1] )
    outdat <- as.list( seq_along( chrlist ) )
    names( outdat ) <- chrlist
    for( chr in chrlist ){
      outdat[[ chr ]] <- indat[ indat[,1] == chr, ncol( indat ) ]
    }
    return( outdat )
  }

  if( class( dat.chip ) == "character" )
    dat.chip <- dataFrameToList( read.table( dat.chip ) )

  if( class( dat.input ) == "character" )
    dat.input <- dataFrameToList( read.table( dat.input ) )
  
  if(is.null(chrlist))chrlist <- as.character(names(dat.chip))

  map <- gc <- nm <- NULL
    ##reading M files
    if(!is.null(mfile)){
      if(!is.null(m.suffix)) {
        for(chr in chrlist){
          fname <- paste(mfile,chr,m.suffix,sep="")
          if(file.exists(fname))
            {
              map <- rbind(map,
                           cbind(
                                 chr,
                                 read.table(fname)[,1:2]
                                 )
                           )
              message(paste("Read mappability file",fname,"for",chr,""))
            }else{
              message(paste(fname,"does not exists;",chr,"passed"))
            }
        }
      }else{
        if(file.exists(mfile)){
          map <- read.table(mfile)
          message(paste("Read mappability file",mfile,"for all chromosomes"))
        }else{
          message(paste(mfile,"does not exists; reading for M passed"))
        }
      }
      map[,1] <- as.character(map[,1])
    }
  
    ##reading GC files
    if(!is.null(gcfile)){
      if(!is.null(gc.suffix)) {
        for(chr in chrlist){
          fname <- paste(gcfile,chr,gc.suffix,sep="")
          if(file.exists(fname))
            {
              gc <- rbind(gc,
                           cbind(
                                 chr,
                                 read.table(fname)[,1:2]
                                 )
                           )
              message(paste("Read GC file",fname,"for chromosome",chr,""))
            }else{
              message(paste(fname,"does not exists;",chr,"passed"))
            }
        }
      }else{
        if(file.exists(gcfile)){
          gc <- read.table(gcfile)
          message(paste("Read GC file",gcfile,"for all chromosomes"))
        }else{
          message(paste(gcfile,"does not exists; GC passed"))
        }
      }
      gc[,1] <- as.character(gc[,1])
    }


    ##reading N files
    if(!is.null(nfile)){
      if(!is.null(n.suffix)) {
        for(chr in chrlist){
          fname <- paste(nfile,chr,n.suffix,sep="")
          if(file.exists(fname))
            {
              nm <- rbind(nm,
                           cbind(
                                 chr,
                                 read.table(fname)[,1:2]
                                 )
                           )
              message(paste("Read N file",fname,"for chromosome",chr,""))
            }else{
              message(paste(fname,"does not exists;",chr,"passed"))
            }
        }
      }else{
        if(file.exists(nfile)){
          nm <- read.table(nfile)
          message(paste("Read N file",nfile,"for all chromosomes"))
        }else{
          message(paste(nfile,"does not exists; N passed"))
        }
      }
      nm[,1] <- as.character(nm[,1])
    }

    coord <- Y <- X <- M <- GC <- chrID <- NULL
    message(paste("chromosomes in chiptag contains",paste(chrlist,collapse=",")))
    message(paste("chromosomes in inputtag contains",paste(as.character(names(dat.input)),collapse=",")))
    message(paste("chromosomes in mappability contains",paste(unique(map[,1]),collapse=",")))
    message(paste("chromosomes in GC contains",paste(unique(gc[,1]),collapse=",")))
    message(paste("chromosomes in N contains",paste(unique(nm[,1]),collapse=",")))
    
    for(chr in chrlist){
      chrM <- chrGC <- chrX <- chrN <- chrPos <- NULL
      n <- length(dat.chip[[chr]])
      message(paste("chiptag contains",n,"bins for chromosome",chr))
      if(!is.null(map)){
        if(sum(map[,1]==chr)==0)          {
          message(paste("no mappability for",chr) )
        }else{
          chrM <- map[map[,1]==chr,]
          n <- min(c(n,nrow(chrM)))
          chrPos <- chrM[,2]
        }
      }

      if(!is.null(dat.input)){
        if(! chr %in% as.character(names(dat.input))) {
          message(paste("no input reads for",chr))
        }else{
          chrX <- dat.input[[chr]]
          n <- min(c(n,length(chrX)))
        }
      }

      if(!is.null(gc)){
        if(sum(gc[,1]==chr)==0) {
          message(paste("no gc for",chr))
        }else{
          chrGC <- gc[gc[,1]==chr,]
          n <- min(c(n,nrow(chrGC)))
          chrPos <- chrGC[,2]
        }
      }

      nID <- 1:n
      if(!is.null(nm) & sum(nm[,1]==chr)>0 ){
        chrN <- nm[nm[,1]==chr,]
        n <- min(c(n,nrow(chrN)))
        nID <- which(chrN[1:n,3]!=1)
        chrPos <- chrN[,2]
      }

      if(n>0){
        Y <- c(Y,dat.chip[[chr]][nID])
        if(!is.null(chrPos))          {
          coord <- c(coord,chrPos[nID])
        }else{
          coord <- c(coord,rep(NA,length(nID)))
        }
        if(!is.null(chrM)) {
          M <- c(M,chrM[nID,3])
        }else{
          M <- c(M,rep(NA,length(nID)))
        }
        if(!is.null(chrGC)){
          GC <- c(GC,chrGC[nID,3])
        }else{
          GC <- c(GC,rep(NA,length(nID)))
        }
        if(!is.null(chrX)) {
          X <- c(X,chrX[nID])
        }else{
          X <- c(X,rep(NA,length(nID)))
        }
        chrID <- c(chrID,rep(chr,length(nID)))
      }
      message(paste("Read in values for",n,"bins for chromosome",chr))

    }
 
    newdat=new( "BinData", chrID=chrID,coord=coord,dataType=dataType )
    if(!is.null(Y)) newdat@tagCount=Y
    if(length(na.omit(M))>0) newdat@mappability=M
    if(length(na.omit(GC))>0) newdat@gcContent=GC
    if(length(na.omit(X))>0) newdat@input=X
    return(newdat)
  }
