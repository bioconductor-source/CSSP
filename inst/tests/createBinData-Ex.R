data(tagdat_chip)
data(tagdat_input)
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
  
  write.table( mapdat[[i]], file = paste("map_chr",i,".txt",sep=""),
              sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table( gcdat[[i]], file = paste("gc_chr",i,".txt",sep=""),
              sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table( ndat[[i]], file = paste("n_chr",i,".txt",sep=""),
              sep = "\t", row.names = FALSE, col.names = FALSE)
}
write.table( allmapdat, file = "allmap.txt" , sep = "\t", row.names = FALSE,
            col.names = FALSE )
write.table( allgcdat,file = "allgc.txt" , sep = "\t", row.names = FALSE,
            col.names = FALSE )
write.table( allndat,file = "alln.txt", sep = "\t", row.names = FALSE,
            col.names = FALSE )

bindata1 <- createBinData( dat_chip, dat_input, mfile = "map_",
                          gcfile = "gc_", nfile = "n_", m.suffix = ".txt",
                          gc.suffix = ".txt", n.suffix = ".txt",
                          chrlist = NULL, dataType = "unique" )
bindata2 <- createBinData( dat_chip, dat_input, mfile = "allmap.txt",
                          gcfile="gc_", nfile = "n_", m.suffix = NULL,
                          gc.suffix = ".txt", n.suffix = ".txt",
                          chrlist = NULL, dataType = "unique" )
bindata3 <- createBinData( dat_chip, dat_input, mfile = "map_",
                          gcfile = "allgc.txt", nfile="n_", m.suffix = ".txt",
                          gc.suffix = NULL, n.suffix = ".txt",
                          chrlist = NULL, dataType = "unique")
bindata4 <- createBinData( dat_chip, dat_input, mfile = "map_",
                          gcfile = "gc_", nfile = "alln.txt", m.suffix = ".txt",
                          gc.suffix = ".txt", n.suffix = NULL,
                          chrlist = NULL, dataType = "unique")

for(i in 1:5){
  for(j in c("map_","gc_","n_")){
    file.remove(paste(j,"chr",i,".txt",sep=""))
  }
}
file.remove("allmap.txt")
file.remove("alln.txt")
file.remove("allgc.txt")
