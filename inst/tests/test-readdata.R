# --------------------------------------------------------------------------------------------------------------------
# author: Chandler Zuo <zuo@stat.wisc.edu>
# ---------------------------------------------------------------------------------------------------------------------

library( testthat )
library( stringr )
library(CSSP)

context( "testing read data" )

data(bindpos)
data(tagdat_chip)
data(tagdat_input)

test_that( "Function tag2bin: converting reads into bin-level counts" ,
{
  ##Arrange

  ## Convert all data into bin data
  bindat_chip<- tag2bin(tagdat_chip,binS=100,fragL=100,prob=1)
  expect_identical(names(bindat_chip),names(tagdat_chip))

  ## Convert 50% of the tags to bin data
  bindat_sub <- tag2bin(tagdat_chip,binS=100,fragL=100,prob=0.5)
  expect_identical(prod(sapply(bindat_chip,sum)>sapply(bindat_sub,sum)),1)

  ## The input data is NULL
  expect_error(bindat <- tag2bin(tagdat1,binS=100,fragL=100),"object 'tagdat1' not found")

})

test_that( "Function createBinData: read in M/GC/N files from either seperate files for each chromosome or a single file for all chromosome" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    mapdat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              M=runif(numBins[i],0.9,1)
                              )
    gcdat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              GC=runif(numBins[i],0.5,1)
                              )
    ndat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              N=rbinom(numBins[i],1,0.01)
                              )
    allmapdat <- rbind(allmapdat,
                       cbind(paste("chr",i,sep=""),mapdat[[i]]))
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    allndat <- rbind(allndat,
                       cbind(paste("chr",i,sep=""),ndat[[i]]))

    write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile="n_", m.suffix = ".txt", gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile="allmap.txt", gcfile="gc_", nfile="n_", m.suffix = NULL, gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata3 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="allgc.txt", nfile="n_", m.suffix = ".txt", gc.suffix = NULL, n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata4 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile="alln.txt", m.suffix = ".txt", gc.suffix = ".txt", n.suffix = NULL,   chrlist = NULL, dataType = "unique")

  expect_identical(bindata1,bindata2)
  expect_identical(bindata1,bindata3)
  expect_identical(bindata1,bindata4)

  for(i in 1:5){
    for(j in c("map_","gc_","n_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
})

test_that( "Function createBinData: read in M/GC/N files, with N file missing or incomplete" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    mapdat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              M=runif(numBins[i],0.9,1)
                              )
    gcdat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              GC=runif(numBins[i],0.5,1)
                              )
    allmapdat <- rbind(allmapdat,
                       cbind(paste("chr",i,sep=""),mapdat[[i]]))
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    if(i<=3){
      ndat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              N=rbinom(numBins[i],1,0.01)
                              )
      allndat <- rbind(allndat,
                       cbind(paste("chr",i,sep=""),ndat[[i]]))
      write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
    write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile=NULL, m.suffix = ".txt", gc.suffix = ".txt", n.suffix = NULL,   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile=NULL, m.suffix = ".txt", gc.suffix = ".txt", n.suffix =".txt",   chrlist = NULL, dataType = "unique")

  expect_identical(bindata1,bindata2)

  for(i in 1:5){
    for(j in c("map_","gc_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  for(i in 1:3){
    for(j in c("n_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
  
})

test_that( "Function createBinData: read in M/GC/N files, with M file missing or incomplete" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    if(i<=3){
      mapdat[[i]] <- data.frame(
                                pos=(0:(numBins[i]-1))*100,
                                M=runif(numBins[i],0.9,1)
                                )
      allmapdat <- rbind(allmapdat,
                         cbind(paste("chr",i,sep=""),mapdat[[i]]))
      write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
    gcdat[[i]] <- data.frame(
                             pos=(0:(numBins[i]-1))*100,
                             GC=runif(numBins[i],0.5,1)
                             )
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    ndat[[i]] <- data.frame(
                            pos=(0:(numBins[i]-1))*100,
                            N=rbinom(numBins[i],1,0.01)
                            )
    allndat <- rbind(allndat,
                     cbind(paste("chr",i,sep=""),ndat[[i]]))
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile=NULL, gcfile="gc_", nfile="n_", m.suffix =NULL, gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile=NULL, gcfile="allgc.txt", nfile="n_", m.suffix = ".txt", gc.suffix = NULL, n.suffix =".txt",   chrlist = NULL, dataType = "unique")
  expect_identical(bindata1,bindata2)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="allgc.txt", nfile="n_", m.suffix = ".txt", gc.suffix = NULL, n.suffix =".txt",   chrlist = NULL, dataType = "unique")

  expect_identical(sum(!is.na(bindata1@mappability[bindata1@chrID=="chr4"]))==0,TRUE)
  expect_identical(sum(!is.na(bindata1@mappability[bindata1@chrID=="chr5"]))==0,TRUE)
  expect_identical(sum(is.na(bindata1@mappability[bindata1@chrID=="chr1"]))==0,TRUE)
  expect_identical(sum(is.na(bindata1@mappability[bindata1@chrID=="chr2"]))==0,TRUE)
  expect_identical(sum(is.na(bindata1@mappability[bindata1@chrID=="chr3"]))==0,TRUE)

  for(i in 1:5){
    for(j in c("n_","gc_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  for(i in 1:3){
    for(j in c("map_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
  
})

test_that( "Function createBinData: read in M/GC/N files, with GC file missing or incomplete" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    if(i<=3){
    gcdat[[i]] <- data.frame(
                             pos=(0:(numBins[i]-1))*100,
                             GC=runif(numBins[i],0.5,1)
                             )
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
      mapdat[[i]] <- data.frame(
                                pos=(0:(numBins[i]-1))*100,
                                M=runif(numBins[i],0.9,1)
                                )
      allmapdat <- rbind(allmapdat,
                         cbind(paste("chr",i,sep=""),mapdat[[i]]))
      write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    ndat[[i]] <- data.frame(
                            pos=(0:(numBins[i]-1))*100,
                            N=rbinom(numBins[i],1,0.01)
                            )
    allndat <- rbind(allndat,
                     cbind(paste("chr",i,sep=""),ndat[[i]]))
    write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile=NULL, nfile="n_", m.suffix = ".txt", gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile="allmap.txt", gcfile=NULL, nfile="alln.txt", m.suffix = NULL, gc.suffix = ".txt", n.suffix =NULL,   chrlist = NULL, dataType = "unique")

expect_identical(bindata1,bindata2)

  bindata1 <- createBinData(dat_chip, dat_input, mfile="allmap.txt", gcfile="gc_", nfile="alln.txt", m.suffix = NULL, gc.suffix = ".txt", n.suffix =NULL,   chrlist = NULL, dataType = "unique")

  expect_identical(sum(!is.na(bindata1@gcContent[bindata1@chrID=="chr4"]))==0,TRUE)
  expect_identical(sum(!is.na(bindata1@gcContent[bindata1@chrID=="chr5"]))==0,TRUE)
  expect_identical(sum(is.na(bindata1@gcContent[bindata1@chrID=="chr1"]))==0,TRUE)
  expect_identical(sum(is.na(bindata1@gcContent[bindata1@chrID=="chr2"]))==0,TRUE)
  expect_identical(sum(is.na(bindata1@gcContent[bindata1@chrID=="chr3"]))==0,TRUE)
  
  for(i in 1:5){
    for(j in c("map_","n_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  for(i in 1:3){
    for(j in c("gc_")){
      file.remove(paste(j,"chr",i,".txt",sep=""))
    }
  }
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
  
})

test_that( "Function createBinData: missing chrs in input/M/GC/N" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)

  dat_input[["chr1"]] <- NULL
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    if(i!=3){
    gcdat[[i]] <- data.frame(
                             pos=(0:(numBins[i]-1))*100,
                             GC=runif(numBins[i],0.5,1)
                             )
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
    if(i!=5){
      mapdat[[i]] <- data.frame(
                                pos=(0:(numBins[i]-1))*100,
                                M=runif(numBins[i],0.9,1)
                                )
      allmapdat <- rbind(allmapdat,
                         cbind(paste("chr",i,sep=""),mapdat[[i]]))
      write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
    if(i!=2){
      ndat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              N=rbinom(numBins[i],1,0.01)
                              )
      allndat <- rbind(allndat,
                       cbind(paste("chr",i,sep=""),ndat[[i]]))
      write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile="n_", m.suffix = ".txt", gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile="allmap.txt", gcfile="allgc.txt", nfile="alln.txt", m.suffix = NULL, gc.suffix = NULL, n.suffix =NULL,   chrlist = NULL, dataType = "unique")

  expect_identical(bindata1,bindata2)
  
  expect_equal(sum(!is.na(bindata1@mappability[bindata1@chrID=="chr5"])),0)
  expect_equal(sum(is.na(bindata1@mappability[bindata1@chrID!="chr5"])),0)
  expect_equal(sum(!is.na(bindata1@gcContent[bindata1@chrID=="chr3"])),0)
  expect_equal(sum(is.na(bindata1@gcContent[bindata1@chrID!="chr3"])),0)
  expect_equal(sum(!is.na(bindata1@input[bindata1@chrID=="chr1"])),0)
  expect_equal(sum(is.na(bindata1@input[bindata1@chrID!="chr1"])),0)

  expect_equal(unique(bindata1@chrID),names(dat_chip))
  
  try({
    for(i in 1:5){
      for(j in c("map_","n_","gc_")){
        file.remove(paste(j,"chr",i,".txt",sep=""))
      }
    }
  },silent=TRUE)
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
  
})

test_that( "Function createBinData: missing chrs in input/M/GC/N" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)

  dat_input[["chr1"]] <- NULL
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    if(i!=3){
    gcdat[[i]] <- data.frame(
                             pos=(0:(numBins[i]-1))*100,
                             GC=runif(numBins[i],0.5,1)
                             )
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
    if(i!=5){
      mapdat[[i]] <- data.frame(
                                pos=(0:(numBins[i]-1))*100,
                                M=runif(numBins[i],0.9,1)
                                )
      allmapdat <- rbind(allmapdat,
                         cbind(paste("chr",i,sep=""),mapdat[[i]]))
      write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
    if(i!=2){
      ndat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              N=rbinom(numBins[i],1,0.01)
                              )
      allndat <- rbind(allndat,
                       cbind(paste("chr",i,sep=""),ndat[[i]]))
      write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
    }
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile="n_", m.suffix = ".txt", gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile="allmap.txt", gcfile="allgc.txt", nfile="alln.txt", m.suffix = NULL, gc.suffix = NULL, n.suffix =NULL,   chrlist = NULL, dataType = "unique")

  expect_identical(bindata1,bindata2)
  
  expect_equal(sum(!is.na(bindata1@mappability[bindata1@chrID=="chr5"])),0)
  expect_equal(sum(is.na(bindata1@mappability[bindata1@chrID!="chr5"])),0)
  expect_equal(sum(!is.na(bindata1@gcContent[bindata1@chrID=="chr3"])),0)
  expect_equal(sum(is.na(bindata1@gcContent[bindata1@chrID!="chr3"])),0)
  expect_equal(sum(!is.na(bindata1@input[bindata1@chrID=="chr1"])),0)
  expect_equal(sum(is.na(bindata1@input[bindata1@chrID!="chr1"])),0)

  expect_equal(unique(bindata1@chrID),names(dat_chip))
  
  try({
    for(i in 1:5){
      for(j in c("map_","n_","gc_")){
        file.remove(paste(j,"chr",i,".txt",sep=""))
      }
    }
  },silent=TRUE)
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
  
})

test_that( "Function createBinData: missing chrs in chip" ,
{
  #Arrange
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)

  dat_chip[["chr5"]] <- NULL
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- NULL
  for(i in 1:5){
    gcdat[[i]] <- data.frame(
                             pos=(0:(numBins[i]-1))*100,
                             GC=runif(numBins[i],0.5,1)
                             )
    allgcdat <- rbind(allgcdat,
                       cbind(paste("chr",i,sep=""),gcdat[[i]]))
    write.table(gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
      mapdat[[i]] <- data.frame(
                                pos=(0:(numBins[i]-1))*100,
                                M=runif(numBins[i],0.9,1)
                                )
      allmapdat <- rbind(allmapdat,
                         cbind(paste("chr",i,sep=""),mapdat[[i]]))
      write.table(mapdat[[i]], file = paste("map_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
      ndat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              N=rbinom(numBins[i],1,0.01)
                              )
      allndat <- rbind(allndat,
                       cbind(paste("chr",i,sep=""),ndat[[i]]))
      write.table(ndat[[i]], file = paste("n_chr",i,".txt",sep=""), sep = "\t",  row.names = FALSE, col.names = FALSE)
  }
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allndat,file=paste("alln.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  
  bindata1 <- createBinData(dat_chip, dat_input, mfile="map_", gcfile="gc_", nfile="n_", m.suffix = ".txt", gc.suffix = ".txt", n.suffix = ".txt",   chrlist = NULL, dataType = "unique")
  bindata2 <- createBinData(dat_chip, dat_input, mfile="allmap.txt", gcfile="allgc.txt", nfile="alln.txt", m.suffix = NULL, gc.suffix = NULL, n.suffix =NULL,   chrlist = NULL, dataType = "unique")

  expect_identical(bindata1,bindata2)
  
   expect_equal(unique(bindata1@chrID),names(dat_chip))
  
  try({
    for(i in 1:5){
      for(j in c("map_","n_","gc_")){
        file.remove(paste(j,"chr",i,".txt",sep=""))
      }
    }
  },silent=TRUE)
  file.remove("allmap.txt")
  file.remove("alln.txt")
  file.remove("allgc.txt")
  
})

test_that("readBinFile - structure tests",
{
  dat_chip <- tag2bin(tagdat_chip,binS=100,fragL=100)
  dat_input <- tag2bin(tagdat_input,binS=100,fragL=100)
  dat_chip[[2]] <- dat_chip[[3]] <- dat_chip[[4]] <- dat_chip[[5]] <- NULL
  dat_input[[2]] <- dat_input[[3]] <- dat_input[[4]] <- dat_input[[5]] <- NULL
  
  numBins <- as.integer(runif(5,190,220))
  mapdat <- gcdat <- ndat <- list(1:5)
  allmapdat <- allgcdat <- allndat <- chipdat <- inputdat <- NULL
  for(i in 1){
    gcdat[[i]] <- data.frame(
                             pos=(0:(numBins[i]-1))*100,
                             GC=runif(numBins[i],0.5,1)
                             )
    allgcdat <- rbind(allgcdat,
                      cbind(paste("chr",i,sep=""),gcdat[[i]]))
    mapdat[[i]] <- data.frame(
                              pos=(0:(numBins[i]-1))*100,
                              M=runif(numBins[i],0.9,1)
                              )
    allmapdat <- rbind(allmapdat,
                       cbind(paste("chr",i,sep=""),mapdat[[i]]))
  }
  
  chipdat <- data.frame(
                        chr=rep("chr1",sum(sapply(dat_chip,length))),
                        pos=(0:(length(dat_input[[1]])-1))*100,
                        count=dat_chip[[1]]
                        )
  inputdat <- data.frame(
                         chr=rep("chr1",sum(sapply(dat_input,length))),
                         pos=(0:(length(dat_input[[1]])-1))*100,
                         count=dat_input[[1]]
                         )
  write.table(allmapdat,file=paste("allmap.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(allgcdat,file=paste("allgc.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(chipdat,file="chip.txt",sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(inputdat,file="input.txt",sep="\t",row.names=FALSE,col.names=FALSE)

  dat1 <- readBinFile(c("M","GC","chip","input"),fileName=c("allmap.txt","allgc.txt","chip.txt","input.txt"))
  dat2 <- createBinData(dat_chip,dat_input,mfile="allmap.txt",gcfile="allgc.txt",nfile=NULL)
  expect_equal(dat2@tagCount,dat1$chip)
  expect_equal(dat2@input,dat1$input)
  expect_equal(dat2@mappability,dat1$M)
  expect_equal(dat2@gcContent,dat1$GC)

  dat3 <- readBinFile(c("GC","chip","input","M"),fileName=c("allgc.txt","chip.txt","input.txt","allmap.txt"))
  expect_equal(dat1,dat3)

  expect_error(readBinFile(c("GC","input","M"),fileName=c("allgc.txt","input.txt","allmap.txt")))
  expect_error(readBinFile(c("GC","chip","input"),fileName=c("allgc.txt","chip.txt","input.txt")))

  file.remove("allmap.txt")
  file.remove("allgc.txt")
  file.remove("chip.txt")
  file.remove("input.txt")

})
