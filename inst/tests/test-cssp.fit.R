# --------------------------------------------------------------------------------------------------------------------
# author: Chandler Zuo <zuo@stat.wisc.edu>
# ---------------------------------------------------------------------------------------------------------------------

library( testthat )
library( stringr )
library(CSSP)

context( "testing cssp.fit and cssp.power" )

data(bin.data)
data(bindata.chr1)
data(tagdat_chip)
data(tagdat_input)

test_that( "Function cssp.fit: check the data.frame and BinData data formats are equivalent" ,
{
  ##Arrange

  write.table(cbind("chr1",bindata.chr1[,c(1,4)]),file="allmap.txt",row.names=FALSE,col.names=FALSE,sep="\t")
  write.table(cbind("chr1",bindata.chr1[,c(1,5)]),file="allgc.txt",row.names=FALSE,col.names=FALSE,sep="\t")

  bindata.1 <- createBinData(tag2bin(tagdat_chip,binS=100,fragL=100),tag2bin(tagdat_input,binS=100,fragL=100),mfile="allmap.txt",gcfile="allgc.txt",nfile=NULL,chrlist="chr1")
  ## Convert all data into bin data

  fit1 <- cssp.fit(bindata.1)
  fit2 <- cssp.fit(bindata.chr1)

  expect_equal(fit1@a,fit2@a)
  expect_equal(fit1@mu.input,fit2@mu.input)

  file.remove("allmap.txt")
  file.remove("allgc.txt")
  
})

test_that( "Function cssp.fit: error messages for BinData " ,
{
  ##Arrange

  bindata1 <- new("BinData",coord=bin.data@coord,gcContent=bin.data@gcContent,tagCount=bin.data@tagCount,input=bin.data@input)
  
  ##Assert
  
  expect_error(cssp.fit(bindata1),"Error: mappability is missing")

    ##Arrange

  bindata1 <- new("BinData",coord=bin.data@coord,mappability=bin.data@mappability,tagCount=bin.data@tagCount,input=bin.data@input)
  
  ##Assert
  
  expect_error(cssp.fit(bindata1,useGrid=1),"Error: useGrid must be logical")
  expect_error(cssp.fit(bindata1,method=1),"Error: method must be either 'mde' or 'gem'")

  ##Arrange

  bindata1 <- new("BinData",coord=bin.data@coord,mappability=bin.data@mappability,gcContent=bin.data@gcContent,input=bin.data@input)
  
  ##Assert
  
  expect_error(cssp.fit(bindata1),"Error: ChIP data is missing")

  ##Arrange

  bindata1 <- new("BinData",coord=bin.data@coord,mappability=bin.data@mappability,gcContent=bin.data@gcContent,tagCount=bin.data@tagCount)
  
  ##Assert
  
  expect_error(cssp.fit(bindata1),"Error: input data is missing")

})

test_that( "Function cssp.fit: error messages for data.frame " ,
{
  ##Arrange

  ##Assert
  expect_error(cssp.fit(bindata.chr1[,-4]),"Error: data.frame must contain chip, input and M columns")
  expect_error(cssp.fit(bindata.chr1[,-3]),"Error: data.frame must contain chip, input and M columns")
  expect_error(cssp.fit(bindata.chr1[,-2]),"Error: data.frame must contain chip, input and M columns")

  ##Arrange
  
  ##Assert

  expect_error(cssp.fit(bindata.chr1,method="aaa"),"Error: method must be either 'mde' or 'gem'")
  expect_error(cssp.fit(bindata.chr1,useGrid=1),"Error: useGrid must be logical")
    
})
