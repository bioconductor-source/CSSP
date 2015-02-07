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
	  if((pos[i]+wsize>=start) && (pos[i]-wsize<=end))
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
	  if((pos2[i]>=start) && (pos1[i]<=end))
	    {
	      REAL(posCount)[i]++;
	    }
	}
      
    }
  UNPROTECT(1);
  return (posCount);
}

SEXP peakcount_uniq(SEXP tagCoord,SEXP peakPos1, SEXP peakPos2, SEXP fragLen)
{
	double *coord=REAL(tagCoord), *pos1=REAL(peakPos1), *pos2=REAL(peakPos2), fragL=REAL(fragLen)[0];
	R_len_t n=length(peakPos1),m=length(tagCoord),i,j,start,end,maxpeaklen;
	SEXP posCount, read1, read2;
	PROTECT(posCount=allocVector(REALSXP,n));
	
	// fixed a bug for non-equal length peaks
	maxpeaklen = pos2[0] - pos1[0];
	for(i = 0; i < n; i ++) {
		if(maxpeaklen < pos2[i] - pos1[i]) {
			maxpeaklen = pos2[i] - pos1[i];
		}
	}

	PROTECT(read1=allocVector(REALSXP,maxpeaklen + fragL + 1 ) );
	PROTECT(read2=allocVector(REALSXP,maxpeaklen + fragL + 1 ) );
	
	for(i=0;i<n;i++)   {
		REAL(posCount)[i]=0;
		for( j = 0; j < pos2[0] - pos1[0] + fragL + 1 ; j ++ ){
			REAL(read1)[j] = 0;
			REAL( read2 )[j] = 0;
		}
		for(j=0;j<m;j++)	{
			if(coord[j]<0)   {
				start=-coord[j]-fragL+1;
			  if(start<0) start=0;
			  end=-coord[j];
			  int id = (int)(start - pos1[i] + fragL);
			  if((pos2[i]>=start) && (pos1[i]<=end) && (REAL(read1)[ id ]==0))	  {
				  REAL(read1)[id]++;
				  REAL(posCount)[i]++;
			  }
			} else  {
				start=coord[j];
				end=coord[j]+fragL-1;
				int id = (int)(start - pos1[i] + fragL);
				if((pos2[i]>=start) && (pos1[i]<=end) && (REAL(read2)[ id ]==0) )	  {
					REAL(read2)[id]++;
					REAL(posCount)[i]++;
				}
			}
		}
	}
	UNPROTECT(3);
	return (posCount);
}
