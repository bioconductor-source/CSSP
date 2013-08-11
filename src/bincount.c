#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>

SEXP bindcount_c(SEXP tagCoord,SEXP bindPos, SEXP fragLen, SEXP whs)
{
  double *coord=REAL(tagCoord), *pos=REAL(bindPos), fragL=REAL(fragLen)[0], wsize=REAL(whs)[0];
  R_len_t n=length(bindPos),m=length(tagCoord),i,j,start,end;
  SEXP posCount;
  PROTECT(posCount=allocVector(REALSXP,n));

  for(i=0;i<n;i++)
    {
       REAL(posCount)[i]=0;
      for(j=0;j<m;j++)
	{
	  if(coord[j]<0)
	    {
	      start=-coord[j]-fragL+1;
	      if(start<0) start=0;
	      end=-coord[j];
	    }else{
	    start=coord[j];
	    end=coord[j]+fragL-1;
	  }
	  if(pos[i]+wsize>=start && pos[i]-wsize<=end)
	    {
	      REAL(posCount)[i]++;
	    }
	}
      
    }
  UNPROTECT(1);
  return (posCount);
}

SEXP peakcount_c(SEXP tagCoord,SEXP peakPos1, SEXP peakPos2, SEXP fragLen)
{
  double *coord=REAL(tagCoord), *pos1=REAL(peakPos1), *pos2=REAL(peakPos2), fragL=REAL(fragLen)[0];
  R_len_t n=length(peakPos1),m=length(tagCoord),i,j,start,end;
  SEXP posCount;
  PROTECT(posCount=allocVector(REALSXP,n));

  for(i=0;i<n;i++)
    {
      REAL(posCount)[i]=0;
      for(j=0;j<m;j++)
	{
	  if(coord[j]<0)
	    {
	      start=-coord[j]-fragL+1;
	      if(start<0) start=0;
	      end=-coord[j];
	    }else{
	    start=coord[j];
	    end=coord[j]+fragL-1;
	  }
	  if(pos2[i]>=start && pos1[i]<=end)
	    {
	      REAL(posCount)[i]++;
	    }
	}
      
    }
  UNPROTECT(1);
  return (posCount);
}
