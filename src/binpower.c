#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>

//#define compute_p

//double compute_p(double, int, double);

double compute_p (double *pval, R_len_t m, double q)
{
  R_len_t i,j,k;
  double tmp;
      k=0;
  for(i=1;i<m;i++)
    {
         k++;
      for(j=0;j<i;j++)
	{
	  if(pval[i]<pval[j])
	    {
	      tmp=*(pval+j);
	      *(pval+j)=*(pval+i);
	      *(pval+i)=tmp;
	    }
	}
          if(k>1000)
           {
             k-=1000;
             Rprintf("%d\n",i);
           }
    }
  tmp=0;
  for(i=m-1;i>-1;i--)
    {
//      Rprintf("ordered_pval_%lf\n",pval[i]);
      if(pval[i]/i*m<q) 
	{
	  tmp=pval[i];
	  break;
	}
    }
  return (tmp);
}


SEXP binpower (SEXP size_back, SEXP scale_back, SEXP size_sig, SEXP scale_sig, SEXP prob_bind, SEXP prob_sig, SEXP ite, SEXP fold, SEXP min_count, SEXP beta, SEXP mu_chip)
{
  R_len_t n= length(mu_chip),i,j,k,nite=REAL(ite)[0],ncomp=length(size_sig)/length(mu_chip);
  double *rsize_back=REAL(size_back), *rscale_back=REAL(scale_back), *rscale_sig=REAL(scale_sig), *rsize_sig=REAL(size_sig), *rprob_bind=REAL(prob_bind), *rprob_sig=REAL(prob_sig);
  double rbeta=REAL(beta)[0], *rmu_chip=REAL(mu_chip), rfold=REAL(fold)[0], rmin_count=REAL(min_count)[0];
  double s,xx;
  
  //compute FDR a
  
  double sim_bin[n],tmpcount,tmp,tmprate,p_thre=0,ptail_sig[ncomp],rej_back,binpow;
  SEXP wei_pow;
  PROTECT(wei_pow=allocVector(REALSXP,1));
  
  GetRNGstate();
  Rprintf("Computing FDR=");
  
  for(i=0;i<nite;i++)    {
    k=0;
    for(j=0;j<n;j++)	{
      k++;
      tmp=rbinom(1,rprob_bind[j]);
      if(tmp==1)   {
	tmp=0;
	s=0;
	xx=runif(0,1);
	while(xx>=s && tmp<ncomp){
	  s+=rprob_sig[j+(int)tmp*n];
	  tmp++;
	}
	tmpcount=rnbinom(rsize_sig[j+((int)tmp-1)*n],1/(1+rscale_sig[j+((int)tmp-1)*n]));
      }else{
	tmpcount=rnbinom(rsize_back[j],1/(1+rscale_back[j]));
      }
      sim_bin[j]=pnbinom(tmpcount,rbeta,rbeta/(rbeta+rmu_chip[j]),0,0);
      if(k>1000)	{
	k-=1000;
	Rprintf("simulated 1000 rounds, %3.0f\n",tmpcount);
      }
    }
    p_thre+=compute_p(sim_bin,n,0.05);
  }
  p_thre/=nite;
  //  Rprintf("%lf\n",p_thre);
  
//  Rprintf("Simulating Bin Test\n");  
  
  double w_bin,sum_w=0,sum_pow=0;
  for(i=0;i<n;i++)
    {
      binpow=0;
      for(j=0;j<ncomp;j++){
	ptail_sig[j]=pgamma(fmax2(rmin_count,rfold*rmu_chip[i]),rsize_sig[i+j*n],rscale_sig[i+j*n],0,0);
      }
      rej_back=qnbinom(p_thre,rbeta,rbeta/(rmu_chip[i]+rbeta),0,0);
      for(j=0;j<nite;j++)	{
	tmp=0;
	xx=runif(0,1);
	s=0;
	while(s<=xx && tmp<ncomp){
	  s+=rprob_sig[i+(int)tmp*n];
	  tmp++;
	}
	xx=runif(0,1);
	tmprate=qgamma(ptail_sig[(int)tmp-1]*xx,rsize_sig[i+((int)tmp-1)*n],rscale_sig[i+((int)tmp-1)*n],0,0);
	tmpcount=rpois(tmprate);
	if(tmpcount>=rej_back)	  {
	  binpow++;
	}
	//          Rprintf("simulated_test:%lf\t%lf\n",tmpcount,rej_back);
      }
      binpow/=nite;
      w_bin=0;
      for(j=0;j<ncomp;j++){
	w_bin+=rprob_sig[i+j*n]*ptail_sig[j]*rprob_bind[i];
      }
      sum_w+=w_bin;
      sum_pow+=w_bin*binpow;
    }
  PutRNGstate();
  UNPROTECT(1);
  REAL(wei_pow)[0]=sum_pow/sum_w;
  
  return (wei_pow);
  
}

SEXP binpower_pval (SEXP size_back, SEXP scale_back, SEXP size_sig,SEXP scale_sig, SEXP prob_bind, SEXP prob_sig,SEXP ite, SEXP fold, SEXP min_count, SEXP beta, SEXP mu_chip, SEXP p_val, SEXP prob_zero)
{
  R_len_t n= length(mu_chip),i,j,isig,nite=REAL(ite)[0],ncomp=length(size_sig)/length(mu_chip);
  double *rscale_sig=REAL(scale_sig), *rsize_sig=REAL(size_sig), *rprob_bind=REAL(prob_bind), *rprob_sig=REAL(prob_sig), *rprob_zero=REAL(prob_zero);
  //  double *rsize_back=REAL(size_back), *rscale_back=REAL(scale_back);
  double rbeta=REAL(beta)[0], *rmu_chip=REAL(mu_chip), rfold=REAL(fold)[0], rmin_count=REAL(min_count)[0], p_thre=REAL(p_val)[0];
  double s,xx;
  //compute FDR a
  
  double tmpcount,tmprate,ptail_sig[ncomp],rej_back,binpow;
  SEXP wei_pow,p_tail;
  PROTECT(wei_pow=allocVector(REALSXP,1));
  PROTECT(p_tail=allocVector(REALSXP,n));
  
  GetRNGstate();
//  Rprintf("Computing FDR=");
  
//  Rprintf("Simulating Bin Test\n");  
  
  double w_bin,sum_w=0,sum_pow=0;
  for(i=0;i<n;i++)    {
    binpow=0;
    for(j=0;j<ncomp;j++){
      ptail_sig[j]=pgamma(fmax2(rmin_count,rfold*rmu_chip[i]),rsize_sig[i+j*n],rscale_sig[i+j*n],0,0);
    }
    rej_back=qnbinom(p_thre / ( 1 - rprob_zero[i] ), rbeta, rbeta/(rmu_chip[i]+rbeta),0,0);
    for(j=0;j<nite;j++)	{
      xx=runif(0,1);
      s=0;
      isig=0;
      while(s<=xx && isig<ncomp){
	s+=rprob_sig[i+isig*n];
	isig++;
      }
      xx=runif(0,1);
      tmprate=qgamma(ptail_sig[isig-1]*xx,rsize_sig[i+(isig-1)*n],rscale_sig[i+(isig-1)*n],0,0);
      tmpcount=rpois(tmprate);
      if(tmpcount>=rej_back && tmpcount>=fmax2(rmin_count,rfold*rmu_chip[i]))	    {
	binpow++;
      }
      //          Rprintf("simulated_test:%lf\t%lf\n",tmpcount,rej_back);
    }
    binpow/=nite;
    w_bin=0;
    for(j=0;j<ncomp;j++){
      //      Rprintf("  %f\t%f\t%f\n",rprob_sig[i+j*n],ptail_sig[j],rprob_bind[i]);
      w_bin+=rprob_sig[i+j*n]*ptail_sig[j]*rprob_bind[i];
    }
    sum_w+=w_bin;
    sum_pow+=w_bin*binpow;
    REAL(p_tail)[i]=w_bin;
  }
  PutRNGstate();
  UNPROTECT(2);
  REAL(wei_pow)[0]=sum_pow;
  
  return (wei_pow);
  
}
