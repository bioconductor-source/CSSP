#'@name cssp.fit
#'@title Fit the CSSP Model.
#'
#'@param dat A \link{data.frame} or \link{BinData-class} object containing bin-level chip, input, M and GC information. For the data.frame object, the columns must contain "chip", "input", "M". For BinData object, the slots must contain "tagCount", "input", "M". If "GC" is not provided, model will be fitted without using gc-Content scores.
#'@param method A \link{character} indicating the method of fitting algorithm to be used. "mde" (Default) - minimum distance estimation; "gem" - the generalized EM method.
#'@param p1 The \link{numeric} value for the lower bound for the p-value region where the p-values are assumed to be uniformly distributed. Default: 0.5.
#'@param p2 The \link{numeric} value for the upper bound for the p-value region where the p-values are assumed to be uniformly distributed. Default: 0.99.
#'@param beta.init The \link{numeric} value for the initializing the size parameter for the background model of the ChIP sample. If "NULL", the size parameter of the fitted input sample model is used.
#'@param e0.init The \link{numeric} value for initializing parameter e0. Default: 0.9.
#'@param ngc An \link{integer} value for the number of knots used in the spline model for the gc covariate. Default: 9.
#'@param nite An \link{integer} value for the maximum number of iterations taken. Default: 50.
#'@param tol A \link{numeric} value for the tolerance for convergence. Default: 1e-3.
#'@param useGrid The \link{logical} value indicating whether the gridding method is used. If TRUE, the covariate space is grided adaptively. This trims down the sample size for fitting the regression model when the data contains too many observations, and is suggested for genome-wide analysis. Default: FALSE.
#'@param nsize A \link{numeric} value for the number of bins to be randomly chosen in estimating the normalizatiing parameters. If Null (default), all bins are used in normalization. For genome wide analysis, nsize=5000 is suggested.
#'@param ncomp A \link{numeric} value for the number of signal components.
#'@return \link{CSSPFit-class} A CSSPFit object.
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( bin.data )
#'cssp.fit( bin.data )
#'cssp.fit( bin.data, method = "gem" )
#'data( bindata.chr1 )
#'cssp.fit( bindata.chr1 )
#'cssp.fit( bindata.chr1, method = "gem", ngc = 1 )
#'@name cssp.fit
#'@aliases cssp.fit,data.frame-method,BinData-method
#'@docType methods
#'@rdname cssp.fit-methods
#'@export
setGeneric("cssp.fit",
           function(dat,...)
           standardGeneric("cssp.fit")
          )

#'@useDynLib CSSP
#'@rdname cssp.fit-methods
#'@aliases cssp.fit,data.frame-method
setMethod("cssp.fit",
          signature="data.frame",
          definition=function(dat,method="mde",p1=0.5,p2=0.99,beta.init=NULL,e0.init=0.90,ngc=9,nite=50,tol=1e-03,useGrid=FALSE,nsize=NULL,ncomp=2)
          {
            if(prod(c("chip","input","M")%in%names(dat))==0)
              stop("Error: data.frame must contain chip, input and M columns")
            if(!method%in%c("mde","gem"))
              stop("Error: method must be either 'mde' or 'gem'")
            if(!is.logical(useGrid))
              stop("Error: useGrid must be logical")
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
            return(.csspfit(sam.chip,sam.input,map,gc,method,p1,p2,beta.init,e0.init,ngc,nite,tol,m,map.id,useGrid,nsize,ncomp))
          }
          )

#'@useDynLib CSSP
#'@rdname cssp.fit-methods
#'@aliases cssp.fit,BinData-method
setMethod("cssp.fit",
          signature="BinData",
          definition=function(dat,method="mde",p1=0.5,p2=0.99,beta.init=NULL,e0.init=0.90,ngc=9,nite=50,tol=1e-03,useGrid=FALSE,nsize=NULL,ncomp=2)
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
            m <- length(dat@tagCount)
            map.id <- which(dat@mappability>0)
            sam.chip <- as.numeric(dat@tagCount[map.id])
            sam.input <- as.numeric(dat@input[map.id])
            gc <- NULL
            if(length(dat@gcContent)>0) gc <- dat@gcContent[map.id]
            map <- dat@mappability[map.id]
            return(.csspfit(sam.chip,sam.input,map,gc,method,p1,p2,beta.init,e0.init,ngc,nite,tol,m,map.id,useGrid,nsize,ncomp))
          }
          )

.csspfit <- function(sam.chip,sam.input,map,gc,method,p1,p2,beta.init,e0.init,ngc,nite,tol,m,map.id,useGrid,nsize,ncomp)
  {
    n <- length(sam.chip)

    ## fit the input sample
    if(useGrid)
      {
        input.fit <- .gridFit(sam.input,map,gc,ngc)
      }else{
        input.fit <- .allFit(sam.input,map,gc,ngc)
      }        
    
    mu.input <- input.fit$predict

    t0 <- Sys.time()
    message("===Computing the size parameter for the input sample===")
    
    m1 <- mean(mu.input)
    m2 <- mean((sam.input-mu.input)^2)
    alpha <- m1^2/(m2-m1)
    if(alpha<0) alpha <- 10
    if(alpha==Inf)alpha <- mean(sam.input)^2/mean(sam.input^2)
    
    if(method=="mde")
      {
        stepsizealpha <- alpha/10
        drctalpha <- 0
        n1alpha <- n
        n2alpha <- 0
        if(is.null(nsize))
          {
            nsize <- n
            sub.ind <- 1:length(sam.input)
          }else{
            sub.ind <- sample(1:length(sam.input),nsize,replace=FALSE)
          }
        sub.input <- sam.input[sub.ind]
        sub.mu <- mu.input[sub.ind]
        sub.chip <- sam.chip[sub.ind]
        
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
            pval1 <- pnbinom(sub.input,size=alpha,mu=sub.mu,lower.tail=FALSE)
            pval2 <- pnbinom(sub.input-1,size=alpha,mu=sub.mu,lower.tail=FALSE)
            pval.null <- pval2+(pval1-pval2)*runif(nsize)
            pval <- sort(pval.null)
            rsd <- pval[1:nsize]-pval[1]-(1-pval[1])/(nsize-1)*(0:(nsize-1))
            n1alpha <- sum(rsd[1:as.integer(nsize/2)]<0)
            n2alpha <- sum(rsd[-(1:as.integer(nsize/2))]<0)
            if(length(drctalpha)>=20)break
          }
      }

    message(paste("Time elapsed:",Sys.time()-t0))
    
    ## fit the chip sample background
    ## iterate to find optimal e0 and beta

    
    lambday <- sum(sam.chip)
    lambdax <- sum(sam.input)
    if(!is.null(beta.init) & method=="mde"){beta <- beta.init}else{beta<-alpha}
    e0 <- e0.init
    pi0 <- 0.9
    
    if(method=="mde")
      {
        t0 <- Sys.time()
        message("===Estimating normalizing parameters for the ChIP sample using MDE===")
        stepsize <- (1-e0)/10
        stepsizebeta <- beta/10
        drctbeta <- 0
        drct <- 0
        dis <- c(2,1)
        i <- as.integer(p1*nsize)
        i.up <- as.integer(p2*nsize)
        pi0 <- 0.5
        betalist <- beta
        e0list <- e0
        
        while((length(drctbeta)<nite & length(drct)<nite & abs(dis[length(dis)]-dis[length(dis)-1])>=tol) | pi0>=1)
          {
            mu.chip <- mu.input*lambday*e0/lambdax
            sub.mu <- mu.chip[sub.ind]
            pval1 <- pnbinom(sub.chip,size=beta,mu=sub.mu,lower.tail=FALSE)
            pval2 <- pnbinom(sub.chip-1,size=beta,mu=sub.mu,lower.tail=FALSE)
            pval <- pval2+(pval1-pval2)*runif(length(sub.chip))
            rsd <- sort(pval)[i:i.up]-sort(pval)[i]-(sort(pval)[i.up]-sort(pval)[i])/(i.up-i)*(0:(i.up-i))
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
            stepsize <- min(c(stepsize,(1-e0)/10))
            mu.chip <- mu.input*lambday*e0/lambdax
            sub.mu <- mu.chip[sub.ind]
            pval1 <- pnbinom(sub.chip,size=beta,mu=sub.mu,lower.tail=FALSE)
            pval2 <- pnbinom(sub.chip-1,size=beta,mu=sub.mu,lower.tail=FALSE)
            pval <- pval2+(pval1-pval2)*runif(nsize)
            rsd <- sort(pval)[i:i.up]-sort(pval)[i]-(sort(pval)[i.up]-sort(pval)[i])/(i.up-i)*(0:(i.up-i))
            n1 <- sum(rsd<0)
            n2 <- sum(rsd[1:as.integer(length(rsd)/2)]<0)
            dis <- c(dis,max(abs(rsd)))
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
                if(drct[length(drctbeta)]*drct[length(drctbeta)-1]==-1)
                  {
                    stepsizebeta <- stepsizebeta/2
                  }
              }
            stepsizebeta <- min(c(stepsizebeta,beta/20))
            pi0 <- 1-(i-sort(pval)[i]/(sort(pval)[i.up]-sort(pval)[i])*(i.up-i))/nsize
            betalist <- c(betalist,beta)
            e0list <- c(e0list,e0)
          }
        e0 <- e0list[which.min(dis)-1]
        beta <- betalist[which.min(dis)-1]
        message(paste("Time elapsed:",Sys.time()-t0))
      }
    
    ## estimate the signal distributions as nbinom
    t0 <- Sys.time()
    
    message("===EM algorithm for estimating the signal parameters===")
    
    pi1 <- 1-pi0
    
    mu.chip <- mu.input*e0*lambday/lambdax
    if(method=="gem")     {
      pval <- pnbinom(sam.chip,size=alpha+mu.input,mu=(alpha+sam.input)/(alpha+mu.input)*mu.chip,lower.tail=FALSE)
    }
    resid.sig <- sort((sam.chip-mu.chip)[order(pval)[1:as.integer(n*pi1)]],decreasing=TRUE)
    mean.sig <- rep(mean(resid.sig),ncomp)
    size.sig <- rep(1.2,ncomp)

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
        
        if(prod(apply(em.track[nrow(em.track)-0:1,],2,function(x)abs(x[1]-x[2])/x[1])<=0.05)==1)break
        if(sum(em.track[nrow(em.track),]<0)>0)break
        if(pi0>1 | beta<0 | sum(size.sig<0)>0)stop("error: fitting algorithm is nonconvergent")
        ## E step
        ## use non-convoluted parameter for signal
        prob.back <- dnbinom(sam.chip,size=beta,mu=mu.chip)
        prob.sig <- prob.g <- matrix(0,nrow=n,ncol=ncomp)
        for(i in seq_len(ncomp)){
          prob.sig[,i] <- dnbinom(sam.chip,size=size.sig[i],mu=mean.sig[i])
        }
        for(i in seq_len(ncomp)){
          prob.g[,i] <- p.sig[i]*prob.sig[,i]/as.vector(prob.sig%*%p.sig)
          prob.g[is.na(prob.g[,i]),i] <- p.sig[i]
        }
        
        prob.z <- 1/(1+mean(1-prob.z)/mean(prob.z)*prob.back/as.vector(prob.sig%*%p.sig))
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
          size.sig[i] <- mean.sig[i]/((var.sig[i])/mean.sig[i]-1)
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
    sub.mu <- mu.chip[sub.ind]
    sub.chip <- sam.chip[sub.ind]
    pval1 <- pnbinom(sub.chip,size=beta,mu=sub.mu,lower.tail=FALSE)
    pval2 <- pnbinom(sub.chip-1,size=beta,mu=sub.mu,lower.tail=FALSE)
    pval <- pval2+(pval1-pval2)*runif(length(sub.chip))
    pval <- sort(pval)
    
    message(paste("Time elapsed:",Sys.time()-t0))
    
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
        post.p.sig=prob.g,
        post.p.bind=prob.z,
        post.shape.sig=post.size.sig,
        post.scale.sig=post.scale.sig,
        post.shape.back=post.size.back,
        post.scale.back=post.scale.back,
        n=n,
        k=ncomp,
        map.id=map.id,
        pvalue=pval
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
