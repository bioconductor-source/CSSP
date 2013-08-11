#'@name cssp.sim
#'@title Simulate bin binding intensities according to the posterior distributions of the fitted CSSP model.
#'
#'@param x A \link{numeric} value for the sequencing depth of the ChIP sample at which the new binding intensities at simulated.
#'@param fit A \link{CSSPFit-class} class object describing the CSSP model.
#'@return A \link{list} object containing
#'
#'\tabular{ll}{
#'chip \tab A \link{numeric} vector for the binding intensities for the ChIP sample.\cr
#'
#'bind \tab A \link{numeric} vector for the simulated binding regions.\cr
#'
#'bind.sig \tab A \link{numeric} vector for the signal component for each bin.\cr
#'}
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( sampleFit )
#'cssp.sim( fit = sampleFit, x = sampleFit@@lambday*0.1 )
#'@export
#'@docType methods
#'@rdname cssp.sim-methods
setGeneric("cssp.sim",
           function(fit,...)
           standardGeneric("cssp.sim")
           )

#'@rdname cssp.sim-methods
#'@aliases cssp.sim,CSSPFit-method
setMethod("cssp.sim",
          signature="CSSPFit",
          definition=function(fit,x=fit@lambday)
          {
            n <- fit@n
            pois.mean <- matrix(0,nrow=n,ncol=fit@k)
            for(i in seq_len(fit@k)){
              pois.mean[,i] <- rgamma(n,shape=fit@post.shape.sig[,i],scale=fit@post.scale.sig[,i])*x/fit@lambday
            }
            bind.id <- rbinom(n,size=1,prob=fit@post.p.bind)
            .rmultinom <- function(x){
              y <- x[1]
              x <- x[-1]
              for(i in seq_along(x)){
                if(sum(x[1:i]<=y)) return(i)
              }
              return(length(x))
            }
            bind.sig <- apply(cbind(runif(n),fit@post.p.sig),1,.rmultinom)
            post.mean <- rgamma(n,shape=fit@post.shape.back,scale=fit@post.scale.back)*x/fit@lambday
            for(i in seq_len(fit@k)){
              post.mean[bind.id*bind.sig==i] <- pois.mean[,i][bind.id*bind.sig==i]
            }
            return(list(chip=post.mean,bind=bind.id,bind.sig=bind.sig))
          }
          )

#'@name callpeak
#'@title Call enriched bins based on the CSSP model.
#'
#'@param fit A \link{CSSPFit-class} object containing the fitted CSSP model.
#'@param chip A \link{numeric} vector containing the bin counts for the ChIP sample.
#'@param depth A \link{numeric} value for the sequencing depth corresponding to the ChIP sample of the "chip" argument. If not provided, sequencing depth of "fit" is used.
#'@param fold A \link{numeric} value for the fold change threshold for peak calling.
#'@param min.count A \link{numeric} value for the minimum ChIP count threshold for peak calling.
#'@param qval A \link{numeric} value for the false-discovery rate to be controlled. Default: 0.05.
#'@param method A \link{character} value. By default, "min.count" is used to threshold the ChIP bin counts. If 'method=="post"', "min.count" is used to threshold the posterior bin-level poisson intensities.
#'@return A \link{numeric} vector of locations for binding bins. 
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( sampleFit )
#'data( bin.data )
#'callpeak( sampleFit, chip = bin.data@@tagCount, fold = 1, min.count = 0 )
#'@export
#'@docType methods
#'@rdname callpeak-methods
setGeneric("callpeak",
           function(fit,...)
           standardGeneric("callpeak")
          )

#'@rdname callpeak-methods
#'@aliases callpeak,CSSPFit-method
setMethod("callpeak",
          signature="CSSPFit",
          definition=function(fit,chip,fold=1.8,min.count=0,qval=0.05,method="",depth=fit@lambday)
          {
            pval <- pnbinom(chip[fit@map.id],mu=fit@mu.chip*depth/fit@lambday,size=fit@b,lower=F)
            thr <- apply(cbind(fold*fit@mu.chip*depth/fit@lambday,min.count),1,max)
            if(method=="post")
              {
                p <- matrix(0,nrow=fit@n,ncol=fit@k)
                for(i in seq_len(fit@k)){
                  p[,i] <- pgamma(thr,shape=fit@post.shape.sig[,i],scale=fit@post.scale.sig[,i]*depth/fit@lambday,lower=F)
                }
#                p0 <- pgamma(thr,shape=fit@post.shape.back,scale=fit@post.scale.back*depth/fit@lambday,lower=F)
                                        #            p.thr <- fit@post.p.bind*(fit@post.p.sig1*p1+fit@post.p.sig2*p2)+(1-fit@post.p.bind)*p0
                p.thr <- fit@post.p.sig%*%p
                return(fit@map.id[which(pval*fit@n/rank(pval)<=qval & p.thr>0.5)])
              }else{
                return(fit@map.id[which(pval*fit@n/rank(pval)<=qval & chip[fit@map.id]>=thr)])
              }
          }
         )

#'@name fit.freq
#'@title Compute the estimated frequency for ChIP counts based on the CSSP model.
#'
#'@param fit A \link{CSSPFit-class} object for the fitted CSSP model.
#'@param chip A \link{numeric} vector of ChIP sample bin counts.
#'@return A \link{data.frame} object containing
#' \tabular{l}{
#' count The counts of each bin.\\
#' freq The ChIP data frequency at this count value.\\
#' freq.est  The estimated frequency using the posterior distributions of the bin-level poisson intensities.}
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( sampleFit )
#'data( bin.data )
#'fit.freq( sampleFit, chip = bin.data@@tagCount )
#'@export
#'@docType methods
#'@rdname fit.freq-methods
setGeneric("fit.freq",
           function(fit,...)
           standardGeneric("fit.freq")
          )

#'@rdname fit.freq-methods
#'@aliases fit.freq,CSSPFit-method
setMethod("fit.freq",
          signature="CSSPFit",
          definition=function(fit,chip)
          {
            if(length(chip)<fit@n)
              stop("chip must have length at least fit@n")
            freq.chip <- matrix( 0, nrow = as.integer( quantile( chip, 0.99 ) ), ncol = 3 )
            for(i in ( seq_len(quantile(chip,0.99) )-1) )
              {
                post.dnbinom <- rep(0,fit@n)
                for(j in seq_len(fit@k)){
                  post.dnbinom <- post.dnbinom +fit@post.p.sig[,j]*dnbinom(i,size=fit@post.shape.sig[,j],mu=fit@post.scale.sig[,j]*fit@post.shape.sig[,j])
                }
                freq.chip[i,] <- c(i,
                                   sum(chip[fit@map.id]==i),
                                   sum(
                                       dnbinom(i,size=fit@post.shape.back,mu=fit@post.shape.back*fit@post.scale.back)*(1-fit@post.p.bind)+fit@post.p.bind*post.dnbinom  )
                                   )
              }
            freq.chip <- data.frame(freq.chip)
            names(freq.chip) <- c("count","freq","freq.est")
            return(freq.chip)
          })
            
#'@name qBBT
#'@title Compute the quantile estimate for the bin-level poisson parameters.
#'
#'@param fit A \link{CSSPFit-class} object for the CSSP model.
#'@param prob A \link{numeric} value for the percentile level of bin-level poisson parameters.
#'@param depth A \link{numeric} value for the sequencing depth at which the quantile is evaluated.
#'@param lower A \link{logical} value. If TRUE, the lower quantile is computed. If FALSE (Default), the upper quantile is computed.
#'@return A \link{numeric} value for the percentile of bin-level poisson parameters. 
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( sampleFit )
#'qBBT( sampleFit, prob = 0.99, depth = sampleFit@@lambday*0.1 )
#'@export
#'@docType methods
#'@rdname qBBT-methods
setGeneric("qBBT",
           function(fit,...)
           standardGeneric("qBBT")
          )

#'@rdname qBBT-methods
#'@aliases qBBT,CSSPFit-method
setMethod("qBBT",
          signature="CSSPFit",
          definition=function(fit,prob,depth=fit@lambday,lower=FALSE)
          {
            if(lower) prob <- 1-prob
            ur <- lr <- quantile(fit@post.shape.back*fit@post.scale.back*depth/fit@lambday,prob)
            while(pBBT(fit,ur,depth=depth)<=prob)
              {
                ur <- ur+fit@mean.sig[1]/2
              }
            while(pBBT(fit,lr,depth=depth)>=prob)
              {
                lr <- lr-fit@mean.sig[1]/2
              }
            if(lr<0)return("error: lr < 0")
            return(uniroot(function(x)pBBT(fit,x,depth=depth)-prob,lower=lr,upper=ur)$root)
          }
         )

#'@name pBBT
#'@title Compute the cumulative probability of the bin-level poisson parameters.
#'
#'@param fit A \link{CSSPFit-class} object for the CSSP model.
#'@param x A \link{numeric} value for the percentile level of bin-level poisson parameters.
#'@param depth A \link{numeric} value for the sequencing depth at which the probability is estimated.
#'@param lower A \link{logical} value. If TRUE, the lower quantile is computed. If FALSE (Default), the upper quantile is computed.
#'@return A \link{numeric} value for the cumulative distribution of bin-level poisson parameters. 
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( sampleFit )
#'pBBT( sampleFit, x = 10 )
#'@export
#'@docType methods
#'@rdname pBBT-methods
setGeneric("pBBT",
           function(fit,...)
           standardGeneric("pBBT")
          )

#'@rdname pBBT-methods
#'@aliases pBBT,CSSPFit-method
setMethod("pBBT",
          signature="CSSPFit",
          definition=function(fit,x,depth=fit@lambday,lower=TRUE)
          {
            post.pgamma <- rep(0,fit@n)
            for(i in seq_len(fit@k)){
              post.pgamma <- post.pgamma + pgamma(x,shape=fit@post.shape.sig[,i],scale=fit@post.scale.sig[,i]*depth/fit@lambday)*fit@post.p.sig[,i]*fit@post.p.bind
            }
            prob <- mean(pgamma(x,shape=fit@post.shape.back,scale=fit@post.scale.back*depth/fit@lambday)*(1-fit@post.p.bind)+post.pgamma)

            if(!lower) prob <- 1-prob
            return(prob)
          }
          )
