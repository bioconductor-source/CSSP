#'@name cssp.fit
#'@title Fit the CSSP Model.
#'
#'@param dat A \link{data.frame} or \link{BinData-class} object containing bin-level chip, input, M and GC information. For the data.frame object, the columns must contain "chip", "input", "M". For BinData object, the slots must contain "tagCount", "input", "M". If "GC" is not provided, model will be fitted without using gc-Content scores.
#'@param method A \link{character} indicating the method of fitting algorithm to be used. "mde" (Default) - minimum distance estimation; "gem" - the generalized EM method.
#'@param p1 The \link{numeric} value for the lower bound for the p-value region where the p-values are assumed to be uniformly distributed. Default: 0.5.
#'@param p2 The \link{numeric} value for the upper bound for the p-value region where the p-values are assumed to be uniformly distributed. Default: 0.99.
#'@param beta.init The \link{numeric} value for the initializing the size parameter for the background model of the ChIP sample. If "NULL", the size parameter of the fitted input sample model is used.
#'@param e0.init The \link{numeric} value for initializing parameter e0. Default: 0.9.
#'@param e0.lb The \link{numeric} value for the lower bound of parameter e0. Default is 0.5. This parameter is recommended to be set according to the p-value plot.
#'@param ngc An \link{integer} value for the number of knots used in the spline model for the gc covariate. Default: 9.
#'@param nite An \link{integer} value for the maximum number of iterations taken. Default: 50.
#'@param tol A \link{numeric} value for the tolerance for convergence. Default: 1e-3.
#'@param useGrid A \link{logical} value indicating whether the gridding method is used. If TRUE, the covariate space is grided adaptively. This trims down the sample size for fitting the regression model when the data contains too many observations, and is suggested for genome-wide analysis. Default: FALSE.
#'@param nsize A \link{numeric} value for the number of bins to be randomly chosen in estimating the normalizatiing parameters. If Null (default), all bins are used in normalization. For genome wide analysis, nsize=5000 is suggested.
#'@param ncomp A \link{numeric} value for the number of signal components.
#'@param nonpa A \link{logical} value indicating whether a nonparametric model for the background ChIP sample and the input sample is fitted.
#'@param zeroinfl A \link{logical} value indicating whether a zero-inflated negative binomial model is fitted for the ChIP background.
#'@param seed A \link{numeric} value for the seed of generating random variables. Default: NULL. Users should specify this value for generating exactly reproducible results.
#'@details The current version of cssp.fit has implemented the following method.\cr
#'The "method" argument specifies the method to estimate the normalization models for the ChIP background from the input data. "mde" uses minimum distance estimation, "gem" uses generalized E-M estimation.\cr
#'The 'nonpa' argument specifies whether a glm model is used. If "nonpa" is FALSE, a GLM is used to fit the input data. If "nonpa" is TRUE, the mean response within each grid is taken as the predict. These two arguments enables the analysis for genome-wide data. In this case, "nsize" grids are used.\cr
#'If "nonpa" is FALSE, then "useGrid" specifies whether the covariate space is grided adaptively, and the mean values within each grid is used for regression.\cr
#'If "nonpa" is TRUE, "zeroinfl" specifies whether a zero-inflation model for the background is used. This is useful for low-depth ChIP data, where too many bins have zero count.
#'@return \link{CSSPFit-class} A CSSPFit object.
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( bin.data )
#'cssp.fit( bin.data )
#'cssp.fit( bin.data, method = "gem" )
#'data( bindata.chr1 )
#'cssp.fit( bindata.chr1 )
#'cssp.fit( bindata.chr1, method = "gem", ngc = 1 )
#'@aliases cssp.fit,data.frame-method,BinData-method
#'@docType methods
#'@rdname cssp.fit-methods
#'@export
setGeneric("cssp.fit",
           function(dat,method="mde",p1=0.5,p2=0.99,beta.init=NULL,e0.init=0.90,e0.lb=0.5,ngc=9,nite=50,tol=0.01,useGrid=FALSE,nsize=NULL,ncomp=2,nonpa=FALSE,zeroinfl=FALSE,seed=NULL)
           standardGeneric("cssp.fit")
          )

#'@useDynLib CSSP
#'@rdname cssp.fit-methods
#'@aliases cssp.fit,data.frame-method
setMethod("cssp.fit",
          signature="data.frame",
          definition=function(dat,method="mde",p1=0.5,p2=0.99,beta.init=NULL,e0.init=0.90,e0.lb=0.5,ngc=9,nite=50,tol=0.01,useGrid=FALSE,nsize=NULL,ncomp=2,nonpa=FALSE,zeroinfl=FALSE,seed=NULL)
          {
            if(prod(c("chip","input","M")%in%names(dat))==0)
              stop("Error: data.frame must contain chip, input and M columns")
            if(!method%in%c("mde","gem"))
              stop("Error: method must be either 'mde' or 'gem'")
            if(!is.logical(useGrid))
              stop("Error: useGrid must be logical")
            if(!is.logical(nonpa))
              stop("Error: nonpa must be logical")
            if(!is.logical(zeroinfl))
              stop("Error: zeroinfl must be logical")
            if( !nonpa & zeroinfl ){
                message( "Warning: zeroinfl is set as FALSE for nonpa=FALSE" )
                zeroinfl <- FALSE
            }
            if( nonpa & !useGrid ){
                message( "Warning: useGrid is set as TRUE for nonpa=TRUE" )
                useGrid <- TRUE
            }
            m <- length(dat$chip)
            map.id <- which(dat$M>0)
            sam.chip <- as.numeric(dat$chip[map.id])
            sam.input <- as.numeric(dat$input[map.id])
            if(!is.null(dat$GC))
              {
                gc <- dat$GC[map.id]
              }else{
                gc <- NULL
              }
            map <- dat$M[map.id]
            return(.csspfit(sam.chip,sam.input,map,gc,method,p1,p2,beta.init,e0.init,e0.lb,ngc,nite,tol,m,map.id,useGrid,nsize,ncomp,nonpa,zeroinfl,seed))
          }
          )

#'@useDynLib CSSP
#'@rdname cssp.fit-methods
#'@aliases cssp.fit,BinData-method
setMethod("cssp.fit",
          signature="BinData",
          definition=function(dat,method="mde",p1=0.5,p2=0.99,beta.init=NULL,e0.init=0.90,e0.lb=0.5,ngc=9,nite=50,tol=0.01,useGrid=FALSE,nsize=NULL,ncomp=2,nonpa=FALSE,zeroinfl=FALSE,seed=NULL)
          {
            if(length(dat@mappability)==0) 
              stop("Error: mappability is missing")
            if(length(dat@input)==0)
              stop("Error: input data is missing")
            if(length(dat@tagCount)==0)
              stop("Error: ChIP data is missing")
            if(!method%in%c("mde","gem"))
              stop("Error: method must be either 'mde' or 'gem'")
            if(!is.logical(useGrid))
              stop("Error: useGrid must be logical")
            if(!is.logical(nonpa))
              stop("Error: nonpa must be logical")
            if(!is.logical(zeroinfl))
              stop("Error: nonpa must be logical")
            if( !nonpa & zeroinfl ){
                message( "Warning: zeroinfl is set as FALSE for nonpa=FALSE" )
                zeroinfl <- FALSE
            }
            if( nonpa & !useGrid ){
                message( "Warning: useGrid is set as TRUE for nonpa=TRUE" )
                useGrid <- TRUE
            }
            m <- length(dat@tagCount)
            map.id <- which(dat@mappability>0)
            sam.chip <- as.numeric(dat@tagCount[map.id])
            sam.input <- as.numeric(dat@input[map.id])
            gc <- NULL
            if(length(dat@gcContent)>0) gc <- dat@gcContent[map.id]
            map <- dat@mappability[map.id]
            return(.csspfit(sam.chip,sam.input,map,gc,method,p1,p2,beta.init,e0.init,e0.lb,ngc,nite,tol,m,map.id,useGrid,nsize,ncomp,nonpa,zeroinfl,seed))
          }
          )

.csspfit <- function(sam.chip,sam.input,map,gc,method,p1,p2,beta.init,e0.init,e0.lb,ngc,nite,tol,m,map.id,useGrid,nsize,ncomp,nonpa,zeroinfl,seed)
  {
      if(is.null(seed))
          seed <- Sys.time()
      set.seed(seed)
    n <- length(sam.chip)
    if( is.null( nsize ) ) nsize <- n

    ## fit the input sample
    if( nonpa ){
      input.fit <- .nonpaFit(sam.input,sam.chip,map,gc,ngc)
    }    else if(useGrid) {
      input.fit <- .gridFit(sam.input,map,gc,ngc)
    }    else{
      input.fit <- .allFit(sam.input,map,gc,ngc)
    }
    
    mu.input <- input.fit$predict

    if( zeroinfl & nonpa ){
      zero.input <- input.fit$predict.zero
      zero.chip <- input.fit$predict.zero.chip
    } else {
      zero.input <- zero.chip <- rep( 0, n )
    }

    t0 <- Sys.time()
    message("===Computing the size parameter for the input sample===")

    if( !nonpa ){
      m1 <- mean(mu.input)
      m2 <- mean((sam.input-mu.input)^2)
      alpha <- m1^2/(m2-m1)
    }else{
      alpha <- input.fit$alpha
    }
    
    if(alpha<0) alpha <- 10
    if(alpha==Inf)alpha <- mean(sam.input)^2/mean(sam.input^2)

    ## truncated continuously corrected p-value
    pval.nb.cont.trunc <- function(x){
      rd <- runif(1)
      if( x[1] == 0 )
        return ( NA )
      return(
             ( rd * pnbinom( x[1]-1, mu= x[3], size=x[4], lower.tail = FALSE )+
             ( 1 - rd ) * pnbinom( x[1], mu=x[3], size=x[4], lower.tail = FALSE ) ) / pnbinom( 0, mu = x[3], size = x[4], lower.tail = FALSE )
             )
    }
    
    ## continuously corrected p-value
    pval.nb.cont <- function(x){
      rd <- runif(1)
      return(
             rd * pnbinom( x[1]-1, mu= x[3], size=x[4], lower.tail = FALSE )+
             ( 1 - rd ) * pnbinom( x[1], mu=x[3], size=x[4], lower.tail = FALSE )
             )
    }

    ## Zero-modified density
    pval.mod.w <- function(x){
      if( x[1] == 0 ){
        a = dnbinom( 0, mu = x[3], size = x[4] )
        if( a < x[2] )
          return(
                 a / ( 1 - a ) * ( 1 - x[2] ) / x[2]
                 )
      }
      return( 1 )
    }
    
    if(method=="mde")
      {
        stepsizealpha <- alpha/10
        drctalpha <- 0
        n1alpha <- n
        n2alpha <- 0
        if(is.null(nsize) | nsize > n)
          {
            nsize <- n
            sub.ind <- 1:length(sam.input)
          }else{
            sub.ind <- sample(1:length(sam.input),nsize,replace=FALSE)
          }
        tmp.data <- cbind( sam.input[sub.ind], zero.input[sub.ind], mu.input[sub.ind], alpha)
        nsize <- length( sub.ind )
        
        while(abs(n1alpha-n2alpha)>sqrt(nsize/4))
          {
            if(length(drctalpha)>2)
              {
                if(drctalpha[length(drctalpha)]!=drctalpha[length(drctalpha)-1])
                  {
                    stepsizealpha <- stepsizealpha/2
                  }
              }
            if(n1alpha>n2alpha)
              {
                if(alpha<stepsizealpha)break
                alpha <- alpha-stepsizealpha
                drctalpha <- c(drctalpha,-1)
              }
            if(n1alpha<n2alpha )
              {
                alpha <- alpha+stepsizealpha
                drctalpha <- c(drctalpha,1)
              }
            tmp.data[,4] <- alpha
            pval.null <- apply( tmp.data, 1, pval.nb.cont )
            pval <- sort(pval.null)
            rsd <- pval[1:nsize]-pval[1]-(1-pval[1])/(nsize-1)*(0:(nsize-1))
            n1alpha <- sum(rsd[1:as.integer(nsize/2)]<0)
            n2alpha <- sum(rsd[-(1:as.integer(nsize/2))]<0)
            if(length(drctalpha)>=50)break
          }
      }

    message("Time elapsed:",Sys.time()-t0)
    
    ## fit the chip sample background
    ## iterate to find optimal e0 and beta
    
    lambday <- sum(sam.chip)
    lambdax <- sum(sam.input)
    if(!is.null(beta.init) & method=="mde"){beta <- beta.init}else{beta<-alpha}
    e0 <- e0.init
    pi0 <- 0.9

    pval.w <- function( x ) return ( 1 )
    if( zeroinfl ){
      pval.w <- pval.mod.w
      ##      pval.nb.cont <- pval.nb.cont.trunc
    }
    
    if(method=="mde")
      {
        t0 <- Sys.time()
        message("===Estimating normalizing parameters for the ChIP sample using MDE===")
        stepsize <- min(c(e0, 1-e0) )/10
        stepsizebeta <- beta/10
        drctbeta <- 0
        drct <- 0
        dis <- c(2,1)
        pi0 <- 0.5
        betalist <- rep(beta,2)
        e0list <- rep( e0, 2 )
        pi0list <- rep( pi0, 2 )
        nsize1 <- length(sub.ind)
        mu.chip <- mu.input*lambday*e0/lambdax
        tmp.dat <- cbind( sam.chip[sub.ind], zero.chip[sub.ind], mu.chip[sub.ind], beta)
        
        while((length(drctbeta)<nite & length(drct)<nite & abs(dis[length(dis)])>=tol ) | pi0>=1 | length( drct ) < 20 )
          {
            if( e0 < e0.lb ) break
            mu.chip <- mu.input*lambday*e0/lambdax
            tmp.dat[,3] <- mu.chip[sub.ind]
            pval <- apply( tmp.dat, 1, pval.nb.cont )
            wei <- apply( tmp.dat, 1, pval.w )
            wei <- wei[order(pval)]
            pval <- sort(pval)
            cum.wei <- cumsum( wei ) / sum( wei )
            i <- min( which( pval > p1 ) )
            i.up <- max( which( pval < p2 ) )
            pi0 <- 1-(pval[i]*cum.wei[i.up]-pval[i.up]*cum.wei[i])/(pval[i]-pval[i.up])
            slope <- pval[i]/(cum.wei[i]-1+pi0)
            rsd <- sort(pval)[i:i.up]- slope * (cum.wei[i:i.up] - 1 + pi0)
            n1 <- sum(rsd<0)
            n2 <- sum(rsd[1:as.integer(length(rsd)/2)]<0)
            if(n2<n1*0.5)
              {
                beta <- beta+stepsizebeta
                drctbeta <- c(drctbeta,1)
              }
            if(n2>n1*0.55)
              {
                beta <- beta-stepsizebeta
                drctbeta <- c(drctbeta,-1)
              }
            if(length(drctbeta)>1)
              {
                if(drctbeta[length(drctbeta)]*drctbeta[length(drctbeta)-1]==-1)
                  {
                    stepsizebeta <- stepsizebeta/2
                  }
              }
            stepsizebeta <- min(c(stepsizebeta,beta/20))
            
            tmp.dat[,4] <- beta
            pval <- apply( tmp.dat, 1, pval.nb.cont )
            wei <- apply( tmp.dat, 1, pval.w )
            wei <- wei[order(pval)]
            pval <- sort(pval)
            cum.wei <- cumsum( wei ) / sum( wei )
            i <- min( which( pval > p1 ) )
            i.up <- max( which( pval < p2 ) )
            pi0 <- 1-(pval[i]*cum.wei[i.up]-pval[i.up]*cum.wei[i])/(pval[i]-pval[i.up])
            slope <- pval[i]/(cum.wei[i]-1+pi0)
            rsd <- sort(pval)[i:i.up]- slope * (cum.wei[i:i.up] - 1 + pi0)
            dis <- c(dis,max(abs(rsd)))
            if(max(rsd)+min(rsd)<0 & pi0<1)
              {
                e0 <- e0+stepsize
                drct <- c(drct,1)
              }else{
                e0 <- e0-stepsize
                drct <- c(drct,-1)
              }
            if(drct[length(drct)]*drct[length(drct)-1]==-1)
              {
                stepsize <- stepsize/2
              }
            stepsize <- min(c(stepsize,(1-e0)/10, e0/10))

            betalist <- c(betalist,beta)
            e0list <- c(e0list,e0)
            pi0list <- c(pi0list,pi0)
          }
        dis <- dis[pi0list<1 & e0list > e0.lb ]
        betalist <- betalist[pi0list<1 & e0list > e0.lb ]
        e0list <- e0list[pi0list<1 & e0list > e0.lb]
        e0 <- e0list[which.min(dis)]
        beta <- betalist[which.min(dis)]
        message(paste("Time elapsed:",Sys.time()-t0))
      }
    
    ## estimate the signal distributions as nbinom
    t0 <- Sys.time()
    
    message("===EM algorithm for estimating the signal parameters===")

    pi1 <- 1-pi0

    if( FALSE ){
      if( zeroinfl ){
        pi1 <- ( 1 - pi0 ) * mean( sam.chip == 0 )
      }
    }
    
    pi0 <- 1-pi1
    
    mu.chip <- mu.input*e0*lambday/lambdax
    if(method=="gem")     {
      pval <- pnbinom(sam.chip,size=alpha+mu.input,mu=(alpha+sam.input)/(alpha+mu.input)*mu.chip,lower.tail=FALSE)
    }
    resid.sig <- sort((sam.chip-mu.chip)[order(pval)[1:as.integer(n*pi1)]],decreasing=TRUE)
    mean.sig <- rep(mean(resid.sig),ncomp)
    size.sig <- rep(1.2,ncomp)


    ## prob.infl: (1-prob.infl)*P(ChIP>0)=P_input(ChIP>0), the probability inflation at 0 to match the probability >0
    prob.infl <- 0
    if( zeroinfl )
      prob.infl <- 1- ( 1 - input.fit$predict.zero.chip ) / ( 1 - dnbinom( 0, mu = mu.chip, size = beta ) )
    prob.infl[prob.infl<0] <- 0
    prob.infl[prob.infl>1] <- 1

    p.sig <- rep(0,ncomp)
    p.sig[1]=0.7
    mean.sig[1] <- mean.sig[1]/2
    if(ncomp>=2){
      for(i in 2:ncomp){
        p.sig[i] <- 0.3/(ncomp-1)
        mean.sig[i] <- mean.sig[i]*5*i
      }
    }
    
    em.track <- matrix(rep(0:1,3*ncomp),nrow=2)
    prob.z <- rep(1-pi0,n)
    
    for(k in seq_len(nite) )
      {
        
        if(prod(apply(em.track[nrow(em.track)-0:1,],2,function(x)abs(x[1]-x[2])/x[1])<=tol)==1)break
        if(sum(em.track[nrow(em.track),]<0)>0)break
        if(pi0>1 | beta<0 | sum(size.sig<0)>0)stop("error: fitting algorithm is nonconvergent")
        ## E step
        ## use non-convoluted parameter for signal
        prob.back <- dnbinom(sam.chip,size=beta,mu=mu.chip)
        prob.sig <- prob.g <- matrix(0,nrow=n,ncol=ncomp)
        for(i in seq_len(ncomp)){
          prob.sig[,i] <- dnbinom(sam.chip,size=size.sig[i],mu=mean.sig[i])
        }
        prob.sig.sum <- as.vector( prob.sig %*% p.sig )
        for(i in 1:ncomp){
          prob.g[,i] <- p.sig[i]*prob.sig[,i]/prob.sig.sum
          prob.g[is.na(prob.g[,i]),i] <- p.sig[i]
        }


        ## prob.zero: posterior probability for the inflated component
        ## prob.z: posterior probability for binding
        prob._ <- prob.infl*mean( 1 - prob.z)
        prob._[sam.chip>0] <- 0
        prob.0 <- prob.back * mean( 1 - prob.z ) * ( 1 - prob.infl )
        prob.j <- mean( prob.z ) * prob.sig.sum
        prob.z <- prob.j / ( prob.j + prob.0 + prob._ )
        prob.zero <- prob._ / (prob.j + prob.0 + prob._ )
        prob.zero[prob.zero>1] <- 1
        prob.z[is.na(prob.z)] <- 1-pi0

        for(i in seq_len(ncomp)){
          p.sig[i] <- weighted.mean(prob.g[,i],prob.z)
        }
        p.sig <- p.sig/sum(p.sig)
        
        ## M step
        ##Using moment estimation
        sig.m1 <- sig.m2 <- var.sig <- size.sig <- rep(0,ncomp)
        for(i in seq_len(ncomp)){
          sig.m1[i] <- weighted.mean(sam.chip,prob.g[,i]*prob.z)
          sig.m2[i] <- weighted.mean(sam.chip^2,prob.g[,i]*prob.z)
          mean.sig[i] <- sig.m1[i]
          var.sig[i] <- sig.m2[i]-sig.m1[i]^2
          if( var.sig[i] > mean.sig[i] )
            size.sig[i] <- mean.sig[i]/((var.sig[i])/mean.sig[i]-1)
          else
            size.sig[i] <- 1000
        }
        
        mu.chip <- mu.input*n*weighted.mean(sam.chip[!is.na(prob.z)],1-na.omit(prob.z))/lambdax
        if(method=="gem")
          {
            back.m1 <- weighted.mean(sam.chip,1-prob.z)
            back.m2 <- weighted.mean(sam.chip^2,1-prob.z)
            var.back <- back.m2-back.m1^2
            
            e0 <- back.m1/mean(mu.input)*lambdax/lambday
            
            mu.chip.m1 <- weighted.mean(mu.chip,1-prob.z)
            mu.chip.m2 <- weighted.mean(mu.chip^2,1-prob.z)
            beta <- mu.chip.m2/(back.m2-mu.chip.m2-mu.chip.m1)
            pi0 <- mean(1-prob.z)
            pi1 <- 1-pi0
          }
        em.track <- rbind(em.track,c(p.sig,mean.sig,size.sig))
        if( k %% 10 == 0 )
            message( k, " iterations passed..." )
      }

    post.size.sig <- post.mean.sig <- post.scale.sig <- matrix(0,nrow=n,ncol=ncomp)
    for(i in seq_len(ncomp)){
      post.size.sig[,i] <- size.sig[i]+sam.chip
      post.mean.sig[,i] <- post.size.sig[,i]/(size.sig[i]/mean.sig[i]+1)
      post.scale.sig[,i] <- 1/(size.sig[i]/mean.sig[i]+1)
      post.size.back <- beta+sam.chip
      post.scale.back <- 1/(beta/mu.chip+1)
    }

    sub.ind <- 1:length(sam.input)
    if(!is.null(nsize)){
      sub.ind <- 1:length(sam.input)
      if(nsize<length(sam.input))
        sub.ind <- sample(1:length(sam.input),nsize,replace=FALSE)
    }
    tmp.dat <- cbind( sam.chip[sub.ind], zero.chip[sub.ind], mu.chip[sub.ind], beta)
    pval <- apply( tmp.dat, 1, pval.nb.cont )
    wei <- apply( tmp.dat, 1, pval.w )
    wei <- wei[order(pval)]
    pval <- sort(pval)
    cum.wei <- cumsum( wei ) / sum( wei )
    
    message("Time elapsed:",Sys.time()-t0)

    new("CSSPFit",
        lambdax=lambdax,
        lambday=lambday,
        e0=e0,
        pi0=pi0,
        mu.chip=mu.input*lambday/lambdax*e0,
        mu.input=mu.input,
        a=alpha,
        b=beta,
        mean.sig=mean.sig,
        size.sig=size.sig,
        p.sig=p.sig,
        prob.zero=prob.infl,
        post.p.sig=prob.g,
        post.p.bind=prob.z,
        post.p.zero=prob.zero,
        post.shape.sig=post.size.sig,
        post.scale.sig=post.scale.sig,
        post.shape.back=post.size.back,
        post.scale.back=post.scale.back,
        n=n,
        k=ncomp,
        map.id=map.id,
        pvalue=pval,
        cum.pval=cum.wei
        )
  }

#'@aliases BinData-method
setMethod( "show", "BinData",
          function( object ){
            cat( "class:", class( object ) )
            cat( "length:", length( object@coord ) )
            cat( "chromosomes", sort( unique( as.character( object@chrID ) ) ) )
            cat( "dataType:", object@dataType )
          }
          )

#'@aliases CSSPFit-method
setMethod( "show", "CSSPFit",
          function( object ){
            cat( "class:", class( object ) )
            cat( "length:", object@n )
            cat( "number of components", object@k )
            cat( "e0:", object@e0 )
            cat( "pi0:", object@pi0 )
            cat( showClass( "CSSPFit" ) )
          }
          )

.gridMGC <- function(y,map,gc,ngrid=1000)
  {
    t0 <- Sys.time()
    message("===Gridding the covariate space===")
    res <- .Call("gridMGCmean_c",y,map,gc,ngrid,PACKAGE="CSSP")
    message(paste("===", Sys.time()-t0,"==="))
    t0 <- Sys.time()
    res.mat <- t(matrix(res[2:(4*res[1]+1)],nrow=4))
    return(data.frame(y=res.mat[,1],
                      map=res.mat[,2],
                      gc=res.mat[,3],
                      n=res.mat[,4]))
  }

.gridM <- function(y,map,ngrid=100)
  {
    t0 <- Sys.time()
    message("===Gridding the covariate space===")
    res <- .Call("gridMmean_c",y,map,ngrid,PACKAGE="CSSP")
    message(paste("===", Sys.time()-t0,"===="))
    res.mat <- t(matrix(res[2:(3*res[1]+1)],nrow=3))
    return(data.frame(y=res.mat[,1],
                      map=res.mat[,2],
                      n=res.mat[,3]))
  }

.gridFit <- function(y,map,gc,ngc)
  {
    if(!is.null(gc))
      {
        data.grid <- .gridMGC(y,map,gc)
        if(ngc>0)
          {
            gc_bs <- splines::bs(gc,knots=quantile(gc,prob=seq(0,1,length=ngc+2)[2:(ngc+1)]))
            gc.grid_bs <- splines::bs(data.grid$gc,knots=quantile(data.grid$gc,prob=seq(0,1,length=ngc+2)[2:(ngc+1)]))
          }else{
            gc_bs <- gc
            gc.grid_bs <- data.grid$gc
          }
      }else{
        data.grid <- .gridM(y,map)
      }

    t0 <- Sys.time()
    message("===Constructing the design matrices===")
    
    map_bs <- splines::bs(map,knots=quantile(map[map<0.9],prob=seq(0.1,0.9,0.1)),degree=1)
    map.grid_bs <- splines::bs(data.grid$map,knots=quantile(data.grid$map[data.grid$map<0.9],prob=seq(0.1,0.9,0.1)),degree=1)
    map_high <- (map>0.8)
    map_low <- (map<0.2)
    map.grid_high <- (data.grid$map>0.8)
    map.grid_low <- (data.grid$map<0.2)

    if(!is.null(gc))
      {
        pred.mat <- model.matrix(~map_bs+map_bs*map_high+map_bs*map_low+gc_bs+gc_bs*map_high+gc_bs*map_low)
        fit.mat <- model.matrix(~map.grid_bs+map.grid_bs*map.grid_high+map.grid_bs*map.grid_low+gc.grid_bs+gc.grid_bs*map.grid_high+gc.grid_bs*map.grid_low)
      }else{
        pred.mat <- model.matrix(~map_bs+map_bs*map_high+map_bs*map_low)
        fit.mat <- model.matrix(~map.grid_bs+map.grid_bs*map.grid_high+map.grid_bs*map.grid_low)
      }

    pred.mat <- data.frame(pred.mat)
    pred.mat$y <- NA
    fit.mat <- data.frame(fit.mat)
    fit.mat$y <- data.grid$y
    names(pred.mat) <- names(fit.mat)
    
    message(paste("Time elapsed:",Sys.time()-t0))
    t0 <- Sys.time()
    message("===Fitting glm model based on the mean values of each grid===")

    input.fit <- glm(y~.,data=fit.mat,weights=data.grid$n,family="poisson")
    input.pred <- exp(predict.glm(input.fit,newdata=pred.mat))
    input.pred[input.pred>max(input.fit$fitted.values[input.fit$fitted.values!=Inf])] <- max(input.fit$fitted.values[input.fit$fitted.values!=Inf])
    input.pred[input.pred<min(input.fit$fitted.values[input.fit$fitted.values!=Inf])] <- min(input.fit$fitted.values[input.fit$fitted.values!=Inf])
    message(paste("Time elapsed:",Sys.time()-t0))
    
    return(list(
                fit=input.fit,
                predict=input.pred))
    
  }

.allFit <- function(y,map,gc,ngc)
  {
    t0 <- Sys.time()
    message("===Constructing the design matrices===")
    if(!is.null(gc))
      {
        if(ngc>0)
          {
            gc_bs <- splines::bs(gc,knots=quantile(gc,prob=seq(0,1,length=ngc+2)[2:(ngc+1)]))
          }else{
            gc_bs <- gc
          }
      }
    map_bs <- splines::bs(map,knots=quantile(map[map<0.9],prob=seq(0.1,0.9,0.1)),degree=1)
    map_high <- (map>0.8)
    map_low <- (map<0.2)

    if(!is.null(gc))
      {
        fit.mat <- model.matrix(~map_bs+map_bs*map_high+map_bs*map_low+gc_bs+gc_bs*map_high+gc_bs*map_low)
      }else{
        fit.mat <- model.matrix(~map_bs+map_bs*map_high+map_bs*map_low)
      }

    message(paste("Time elapsed:",Sys.time()-t0))
    t0 <- Sys.time()
    message("===Fitting glm model based on raw values===")

    fit.mat <- data.frame(fit.mat)
    fit.mat$y <- y
    input.fit <- glm(y~.,data=fit.mat,family="poisson")
    input.pred <-input.fit$fitted.values

    input.pred[input.pred>max(input.fit$fitted.values[input.fit$fitted.values!=Inf])] <- max(input.fit$fitted.values[input.fit$fitted.values!=Inf])
    input.pred[input.pred<min(input.fit$fitted.values[input.fit$fitted.values!=Inf])] <- min(input.fit$fitted.values[input.fit$fitted.values!=Inf])
    message(paste("Time elapsed:",Sys.time()-t0))

    return(list(
                fit=input.fit,
                predict=input.pred))
    
  }

.gridMGC0 <- function(x,y,map,gc,ngrid=1000)
  {
    t0 <- Sys.time()
    print("===Gridding the covariate space===")
    res <- .Call("gridMGCnonpa_c",x,y,map,gc,ngrid,PACKAGE="CSSP")
    print(paste("===", Sys.time()-t0,"==="))
    t0 <- Sys.time()
    res.mat <- t(matrix(res[2:(8*res[1]+1)],nrow=8))
    return(list(m=res.mat[,1],
                m2=res.mat[,2],
                m3=res.mat[,3],
                map=res.mat[,4],
                gc=res.mat[,5],
                n0=res.mat[,6],
                n0y=res.mat[,7],
                n=res.mat[,8],
                grid=res[(8*res[1]+2):(8*res[1]+1+length(map))]
                ))
  }

.gridM0 <- function(x,y,map,ngrid=100)
  {
    t0 <- Sys.time()
    print("===Gridding the covariate space===")
    res <- .Call("gridMnonpa_c",x,y,map,ngrid,PACKAGE="CSSP")
    print(paste("===", Sys.time()-t0,"===="))
    res.mat <- t(matrix(res[2:(7*res[1]+1)],nrow=7))
    return(list(m=res.mat[,1],
                m2=res.mat[,2], 
                m3=res.mat[,3], 
                map=res.mat[,4],
                n0=res.mat[,5],
                n0y=res.mat[,6],
                n=res.mat[,7],
                grid=res[(7*res[1]+2):(7*res[1]+1+length(map))]
                ))
  }

.nonpaFit <- function(x,y,map,gc,ngc){
  t0 <- Sys.time()
  message("===Gridding the design space===")
  if(!is.null(gc) )  {
    data.grid <- .gridMGC0(x,y,map,gc)
    if(ngc>0) {
      gc_bs <- splines::bs(gc,knots=quantile(gc,prob=seq(0,1,length=ngc+2)[2:(ngc+1)]))
      gc.grid_bs <- splines::bs(data.grid$gc,knots=quantile(data.grid$gc,prob=seq(0,1,length=ngc+2)[2:(ngc+1)]))
    }else{
      gc_bs <- gc
      gc.grid_bs <- data.grid$gc
    }
  }else{
    data.grid <- .gridM0(x,y,map)
  }
  
  message("Time elapsed:",Sys.time()-t0)

  input.pred <- data.grid$m[data.grid$grid]
  zero.pred <- (data.grid$n0/data.grid$n)[data.grid$grid]
  zero.pred.chip <- (data.grid$n0y/data.grid$n)[data.grid$grid]
  size <- data.grid$m^2/(data.grid$m2-data.grid$m-data.grid$m^2)
  size <- weighted.mean( size[size>0], data.grid$n[size>0] )
  
  input.fit <- list(
                    predict=input.pred,
                    predict.zero = zero.pred,
                    predict.zero.chip = zero.pred.chip,
                    alpha = size)
  return( input.fit )
}
  
