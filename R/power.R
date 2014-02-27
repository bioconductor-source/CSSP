#'@name cssp.power
#'@title Compute the weighted average of bin-wise power conditioning on the fold change and minimal ChIP count requirements.
#'
#'@param x A \link{numeric} value for the sequencing depth of the ChIP sample at which the power is evaluated.
#'@param fit A \link{CSSPFit-class} object for the CSSP model.
#'@param ite A \link{integer} value for the number of iterations used for Monte-Carlo evaluation.
#'@param fold A \link{numeric} value for the fold change threshold.
#'@param min.count A \link{numeric} value for the minimal count threshold.
#'@param useC A \link{logical} value. Whether the function will be evaluated using C. Default: FALSE.
#'@param qval A \link{numeric} value for the q-value for FDR control. Default: 0.05.
#'@return A \link{numeric} value for the weighted average of bin power conditioning on the minimal count and fold change thresholds.
#'@author Chandler Zuo \email{zuo@@stat.wisc.edu}
#'@examples
#'data( sampleFit )
#'cssp.power( sampleFit, x = sampleFit@@lambday*0.1, min.count = 0, fold = 2,
#'useC = TRUE )
#'@export
#'@docType methods
#'@rdname cssp.power-methods
setGeneric("cssp.power",
           function(fit,x,ite=100,fold=1,min.count =10,useC=FALSE,qval=0.05)
           standardGeneric("cssp.power")
           )

#'@useDynLib CSSP
#'@rdname cssp.power-methods
#'@aliases cssp.power,CSSPFit-method
setMethod("cssp.power",
          signature="CSSPFit",
          definition=function(fit,x,ite=100,fold=1,min.count =10,useC=FALSE,qval=0.05)
          {
            n <- length(fit@mu.chip)*ite
            message("===Find the pvalue threshold for FDR control===")
            t0 <- Sys.time()
            a <- rep( 0, ite )
            .rmultinom <- function(y,x){
              for(i in seq_len(length(x))){
                if( sum(x[1:i]) >= y ) return(i)
                return(length(x))
              }
            }
            for(i in seq_len(ite))              {
              pois.mean <- bind.sig <- matrix(0,nrow=n,ncol=fit@k)
              for(j in seq_len(fit@k)){
                pois.mean[,j] <- rgamma(n,shape=fit@post.shape.sig[,j],scale=fit@post.scale.sig[,j])*x/fit@lambday
              }
              bind.id <- rbinom(n,size=1,prob=fit@post.p.bind)
              bind.sig <- apply(cbind(fit@n,fit@post.p.sig),1,function(x).rmultinom(x[1],x[-1]))
              post.mean <- rgamma(n,shape=fit@post.shape.back,scale=fit@post.scale.back)*x/fit@lambday
              infl.id <- rbinom(n,size=1,prob=fit@post.p.zero/(1-fit@post.p.bind))
              post.mean[infl.id] <- 0
              for(j in seq_len(fit@k)){
                post.mean[bind.id==1 & bind.sig==j] <- pois.mean[,j][bind.id==1 & bind.sig==j]
              }
              y.back <- rpois(n,post.mean)
              pval.sim <- pnbinom(y.back,mu=fit@mu.chip*x/fit@lambday,size=fit@b,lower.tail=F) * (1-fit@prob.zero)# + fit@prob.zero
              a[i] <- max(pval.sim[pval.sim/rank(pval.sim)*n<qval] )
            }
            a[a==-Inf] <- 0
            a <- mean(na.omit(a))
            message(paste("Adjusted p-value is",a))
            message(paste("Time elapsed:",Sys.time()-t0))
            t0 <- Sys.time()
            message("===Compute the power===")
            prob.zero <- fit@prob.zero
            if( length( prob.zero ) == 1 )
                prob.zero <- rep( prob.zero, fit@n )
            if(useC==TRUE)              {
              blocksize <- min(c(fit@n/2,10000))
              n.block <- as.integer((fit@n-1)/blocksize)
              bin.pow <- 0
              for(i in seq_len(n.block+1))                {
                if(i<=n.block)                      {
                  id.block <- (1:blocksize)+(i-1)*blocksize
                }else{
                  id.block <- (n.block*blocksize+1):fit@n
                }
                bin.pow <- bin.pow+ .Call("binpower_pval",
                                          fit@post.shape.back[id.block],
                                          fit@post.scale.back[id.block]*x/fit@lambday,
                                          as.vector(fit@post.shape.sig[id.block,]),
                                          as.vector(fit@post.scale.sig[id.block,]*x/fit@lambday),
                                          fit@post.p.bind[id.block],
                                          as.vector(fit@post.p.sig[id.block,]),
                                          ite,
                                          fold,
                                          min.count,
                                          fit@b,
                                          fit@mu.chip[id.block]*x/fit@lambday,
                                          a,
                                          prob.zero[id.block])
              }
              thr <- apply(cbind(fold*fit@mu.chip*x/fit@lambday,min.count),1,max)
              ptail.bind <- matrix(0,nrow=fit@n,ncol=fit@k)
              for(j in seq_len(fit@k)){
                ptail.bind[,j] <- pgamma(thr,shape=fit@post.shape.sig[,j],scale=fit@post.scale.sig[,j]*x/fit@lambday,lower.tail=FALSE)
              }
              w.pow <- sum(as.vector(apply(fit@post.p.sig*ptail.bind,1,sum))*fit@post.p.bind)
              bin.pow <- bin.pow/w.pow
              
            }else{
              crit.val <- qnbinom(1-a/(1-fit@prob.zero),size=fit@b,mu=fit@mu.chip/fit@lambday*x)
              ##   simulate tail counting
              thr <- apply(cbind(fold*fit@mu.chip*x/fit@lambday,min.count),1,max)
              ## the probability for each bin to have strong intensity for each component
              ptail.bind <- matrix(0,nrow=fit@n,ncol=fit@k)
              for(j in seq_len(fit@k)){
                ptail.bind[,j] <- pgamma(thr,shape=fit@post.shape.sig[,j],scale=fit@post.scale.sig[,j]*x/fit@lambday,lower.tail=FALSE)
              }

              ptail.back <- rep(pgamma(thr,shape=fit@post.shape.back,scale=fit@post.scale.back*x/fit@lambday,lower.tail=FALSE),ite)

              ## the strong bin-level intensity for each signal component
              mean.tail.bind <- matrix(0,nrow=fit@n*ite,ncol=fit@k)
              for(j in seq_len(fit@k)){
                mean.tail.bind[,j] <- qgamma(ptail.bind[,j]*runif(fit@n*ite),shape=fit@post.shape.sig[,j],scale=fit@post.scale.sig[,j]*x/fit@lambday,lower.tail=FALSE)
              }
              
              ## choose a signal component for each bin
              pois.tail <- rep(0,fit@n*ite)
              p.sig.mat <- matrix(0,nrow=fit@n*ite,ncol=fit@k)
              for( j in seq_len( fit@k ) ){
                p.sig.mat[,j] <- rep( fit@post.p.sig[,j], ite )
              }
              bind.sig <- apply(cbind(runif(n),p.sig.mat),1,function(x).rmultinom(x[1],x[-1]))
              for(j in seq_len(fit@k)){
                pois.tail[bind.sig==j] <- mean.tail.bind[,j][bind.sig==j]
              }

              ## simulate the bin-level count
              y.tail <- rpois(fit@n*ite,lambda=pois.tail)

              ## the weight
              p.tail <- apply(fit@post.p.sig*ptail.bind,1,sum)*fit@post.p.bind

              ## weight, conditional on the bin is bound
              bin.power <- apply(matrix((y.tail>=crit.val),ncol=ite),1,mean)
              w.power <- rep(0,fit@n*ite)
              for(j in seq_len(fit@k)){
                w.power <- w.power+rep(fit@post.p.sig[,j],ite)*ptail.bind[,j]
              }
              
              bin.pow <- weighted.mean((y.tail>=crit.val&y.tail>=thr)[!is.na(y.tail)],rep(p.tail,ite)[!is.na(y.tail)])
            }
            
            message(paste("Time elapsed:",Sys.time()-t0))
            return(bin.pow)
          }
          )
