#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>

SEXP tag2bin_c(SEXP tagCoord, SEXP fragLen,SEXP binSize,SEXP maxTag)
{
  double *coord=REAL(tagCoord);
  R_len_t fragL=REAL(fragLen)[0], binsize=REAL(binSize)[0], maxpos=REAL(maxTag)[0];
  R_len_t n=length(tagCoord),nbin=maxpos/binsize+1;

  R_len_t startbin,endbin,i,j;

  SEXP binCount;
  PROTECT(binCount=allocVector(REALSXP,nbin));

  for(j=0;j<nbin;j++)
    REAL(binCount)[j]=0;
  
  for(i=0;i<n;i++)
    {
      if(coord[i]<0)
	{
	  startbin=(-coord[i]-fragL+1)/binsize;
	  if(startbin<0) startbin=0;
	  endbin=(-coord[i])/binsize;
	}else{
	startbin=coord[i]/binsize;
	if(startbin<0) startbin=0;
	endbin=(coord[i]+fragL-1)/binsize;
      }
      for(j=startbin;j<endbin+1;j++)
	{
	  REAL(binCount)[j]++;
	}
    }
  UNPROTECT(1);
  
  return (binCount);
}

