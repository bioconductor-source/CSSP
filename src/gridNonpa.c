#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>

SEXP gridMGCnonpa_c(SEXP x,SEXP y, SEXP map, SEXP gc,SEXP ngrid)
{

	double *M=REAL(map),*GC=REAL(gc), *X=REAL(x), *Y=REAL(y);
  R_len_t ng=REAL(ngrid)[0],k,i,j,l,ind[4],n=length(map),newG,counter;

  double Mgrid_l[ng+10],GCgrid_l[ng+10],Ngrid[ng+10],Mgrid_u[ng+10],GCgrid_u[ng+10];
  SEXP res;
  PROTECT(res=allocVector(REALSXP,8*ng+n+50));
  
  Mgrid_l[0]=-0.001;
  Mgrid_u[0]=1;
  GCgrid_l[0]=-0.001;
  GCgrid_u[0]=1;
  k=1;
  Ngrid[0]=n;

  Rprintf("Starting gridding...");
  while(k<ng)
    {
      //find the grid with maximum number of points and which is breakable
      //Rprintf("=====Find the maximum grid=====\n");
      l=0;
      while(((Mgrid_u[l]<=Mgrid_l[l]+0.01 & GCgrid_u[l]>GCgrid_l[l]+0.01) | Mgrid_u[l]+Mgrid_l[l]<0 | GCgrid_u[l]+GCgrid_l[l]<0 )&l<k-1)
	{
	  l++;
	}
      for( i=l+1;i<k;i++)
	{
	  if(Ngrid[i]>Ngrid[l] & (Mgrid_u[i]>Mgrid_l[i]+0.01 | GCgrid_u[i]>GCgrid_l[i]+0.01) & Mgrid_u[l]+Mgrid_l[l]>0 & GCgrid_u[l]+GCgrid_l[l]>0)
	    {
	      l=i;
	    }
	  //Rprintf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,Mgrid_l[i],Mgrid_u[i],GCgrid_l[i],GCgrid_u[i],Ngrid[i]);
	}
      //Rprintf("Maximum grid identified at\t%d\t with count %lf\n",l,Ngrid[l]);
      //Check if breaking this grid is possible
      if(Mgrid_u[l]<=Mgrid_l[l]+0.01 & GCgrid_u[l]<=GCgrid_l[l]+0.01)
	break;
      if(Mgrid_u[l]+Mgrid_l[l]<=0 |GCgrid_u[l]+GCgrid_l[l]<=0)
	break;
      if(Ngrid[l]<=10)
	break;
      //Break this grid into 4 parts and compute the points in each part
      //the topright grid
      Mgrid_u[k]=Mgrid_u[l];
      GCgrid_u[k]=GCgrid_u[l];
      Mgrid_l[k]=(Mgrid_l[l]+Mgrid_u[l])/2;
      GCgrid_l[k]=(GCgrid_l[l]+GCgrid_u[l])/2;
      //the bottomleft grid
      Mgrid_u[l]=Mgrid_l[k];
      GCgrid_u[l]=GCgrid_l[k];
      //the topleft grid
      Mgrid_l[k+1]=Mgrid_l[l];
      Mgrid_u[k+1]=Mgrid_u[l];
      GCgrid_l[k+1]=GCgrid_l[k];
      GCgrid_u[k+1]=GCgrid_u[k];
      //the bottomright grid
      Mgrid_l[k+2]=Mgrid_l[k];
      Mgrid_u[k+2]=Mgrid_u[k];
      GCgrid_l[k+2]=GCgrid_l[l];
      GCgrid_u[k+2]=GCgrid_u[l];

      //compute the number of points in each small grid
      
      //Rprintf("=====New grid counts=====\n");
      for(i=k;i<k+3;i++)
	{
	  Ngrid[i]=0;
	  for(j=0;j<n;j++)
	    {
	      if(M[j]>Mgrid_l[i] & M[j]<=Mgrid_u[i] & GC[j]>GCgrid_l[i] & GC[j]<=GCgrid_u[i])
		{
		  Ngrid[i]=Ngrid[i]+1;
		}
	    }
	  //Rprintf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,Ngrid[i],Mgrid_l[i],Mgrid_u[i],GCgrid_l[i],GCgrid_u[i]);
	}
      Ngrid[l]=Ngrid[l]-Ngrid[k]-Ngrid[k+1]-Ngrid[k+2];
      //Rprintf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",l,Ngrid[l],Mgrid_l[l],Mgrid_u[l],GCgrid_l[l],GCgrid_u[l]);
      
      //delete the 0-count grids
      //number of new grids
      newG=3;
      //if the k to k+2 grids have 0 count
      for(i=1;i<4;i++)
	{
	  if(Ngrid[k+i-1]>0)
	    {
	      ind[i]=k+i-1;
	    }else
	    {
	      ind[i]=-1;
	      newG--;
	    }
	}
      //if the new bottomleft grid has 0 count
      if(Ngrid[l]>0)
	{
	  ind[0]=l;
	}else
	{
	  ind[0]=-1;
	  newG--;
	}
      
      i=0;
      while(ind[i]<0)
	{
	  i++;
	}
      Ngrid[l]=Ngrid[ind[i]];
      Mgrid_l[l]=Mgrid_l[ind[i]];
      Mgrid_u[l]=Mgrid_u[ind[i]];
      GCgrid_l[l]=GCgrid_l[ind[i]];
      GCgrid_u[l]=GCgrid_u[ind[i]];
      //Rprintf("=====New Grids Finally Allocated===\n");
      //Rprintf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",l,Ngrid[l],Mgrid_l[l],Mgrid_u[l],GCgrid_l[l],GCgrid_u[l]);

      j=k;
      while(i<3)
	{
	  i++;
	  if(ind[i]>=0)
	    {
	      Ngrid[j]=Ngrid[ind[i]];
	      Mgrid_l[j]=Mgrid_l[ind[i]];
	      Mgrid_u[j]=Mgrid_u[ind[i]];
	      GCgrid_l[j]=GCgrid_l[ind[i]];
	      GCgrid_u[j]=GCgrid_u[ind[i]];
	      //Rprintf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",j,Ngrid[j],Mgrid_l[j],Mgrid_u[j],GCgrid_l[j],GCgrid_u[j]);
	      j++;
	    }
	}

      k=k+newG;
      //Rprintf("added\t%d\tgrids, new number of grids\t%d\n",newG,k);

    }
  
  Rprintf("finished gridding");

  REAL(res)[0]=k;
  for(i=0;i<k;i++)
    {
      REAL(res)[1+8*i]=0;
      REAL(res)[2+8*i]=0;
      REAL(res)[3+8*i]=0;
      REAL(res)[4+8*i]=0;
      REAL(res)[5+8*i]=0;
      REAL(res)[6+8*i]=0;//number of zero obs of input
      REAL(res)[8+8*i]=Ngrid[i];//number of non-zero obs of input
      REAL(res)[7+8*i]=0;//number of zero obs of chip
      for(j=0;j<n;j++)
	{
	  if(M[j]>Mgrid_l[i] & M[j]<=Mgrid_u[i] & GC[j]>GCgrid_l[i] & GC[j]<=GCgrid_u[i])
	    {
		    REAL(res)[j+8*k+1]=i+1;

		    REAL(res)[1+8*i]+=X[j];
		    REAL(res)[2+8*i]+=X[j]*X[j];
		    REAL(res)[3+8*i]+=X[j]*X[j]*X[j];
		    REAL(res)[4+8*i]+=M[j];
		    REAL(res)[5+8*i]+=GC[j];
		    if( X[j] < 1 )
			    REAL(res)[6+8*i] ++;
		    if( Y[j] < 1 )
			    REAL(res)[7+8*i] ++;
	    }
	}
      //first moment
      REAL(res)[1+8*i] /= Ngrid[i];
      //second moment
      REAL(res)[2+8*i] /= Ngrid[i];
      REAL(res)[3+8*i] /= Ngrid[i];
      REAL(res)[4+8*i] /= Ngrid[i];
      REAL(res)[5+8*i] /= Ngrid[i];
    }

  UNPROTECT(1);
  return(res);

}

SEXP gridMnonpa_c(SEXP x, SEXP y, SEXP map,SEXP ngrid)
{
	double *M=REAL(map), *X=REAL(x), *Y=REAL(y);
	R_len_t ng=REAL(ngrid)[0],k,i,j,l,ind[2],n=length(map),newG,counter;
	
	double Mgrid_l[ng+10],Ngrid[ng+10],Mgrid_u[ng+10];
	SEXP res;
  PROTECT(res=allocVector(REALSXP,7*ng+n+30));
  
  Mgrid_l[0]=-0.001;
  Mgrid_u[0]=1;
  k=1;
  Ngrid[0]=n;

  while(k<ng)
    {
      //find the grid with maximum number of points and which is breakable
      //Rprintf("=====Find the maximum grid=====\n");
      l=0;
      while((Mgrid_u[l]<=Mgrid_l[l]+0.001 | Mgrid_u[l]+Mgrid_l[l]<0)&l<k-1)
	{
	  l++;
	}
      for( i=l+1;i<k;i++)
	{
	  if(Ngrid[i]>Ngrid[l] & Mgrid_u[i]>Mgrid_l[i]+0.001 & Mgrid_u[i]+Mgrid_l[i]>0)
	    {
	      l=i;
	    }
	  //Rprintf("%d\t%lf\t%lf\t%lf\n",i,Mgrid_l[i],Mgrid_u[i],Ngrid[i]);
	}
      //Rprintf("Maximum grid identified at\t%d\t with count %lf\n",l,Ngrid[l]);
      //Check if breaking this grid is possible
      if(Mgrid_u[l]<=Mgrid_l[l]+0.001)
	break;
      if(Mgrid_u[l]+Mgrid_l[l]<=0)
	break;
      if(Ngrid[l]<=10)
	break;
      //Break this grid into 2 intervals
      Mgrid_u[k]=Mgrid_u[l];
      Mgrid_l[k]=(Mgrid_l[l]+Mgrid_u[l])/2;
      //the shrinked old interval
      Mgrid_u[l]=Mgrid_l[k];

      //compute the number of points in each splitted interval
      
      //Rprintf("=====New grid counts=====\n");
      Ngrid[k]=0;
      for(j=0;j<n;j++)
	{
	  if(M[j]>Mgrid_l[i] & M[j]<=Mgrid_u[i])
	    {
	      Ngrid[k]=Ngrid[k]+1;
	    }
	}
      //Rprintf("%d\t%lf\t%lf\t%lf\n",k,Ngrid[k],Mgrid_l[k],Mgrid_u[k]);

      Ngrid[l]=Ngrid[l]-Ngrid[k];
      //Rprintf("%d\t%lf\t%lf\t%lf\n",l,Ngrid[l],Mgrid_l[l],Mgrid_u[l]);
      
      //delete the empty interval
      //number of new intervals
      newG=1;
      //if the k-th interval is empty
      if(Ngrid[k]>0)
	{
	  ind[1]=k;
	}else
	{
	  ind[1]=-1;
	  newG--;
	}
      //if the shrinked old interval is empty
      if(Ngrid[l]>0)
	{
	  ind[0]=l;
	}else
	{
	  ind[0]=-1;
	  newG--;
	}
      
      i=0;
      while(ind[i]<0)
	{
	  i++;
	}
      Ngrid[l]=Ngrid[ind[i]];
      Mgrid_l[l]=Mgrid_l[ind[i]];
      Mgrid_u[l]=Mgrid_u[ind[i]];
      //Rprintf("=====New Grids Finally Allocated===\n");
      //Rprintf("%d\t%lf\t%lf\t%lf\n",l,Ngrid[l],Mgrid_l[l],Mgrid_u[l]);

      j=k;
      if(i<1)
	{
	  i++;
	  if(ind[i]>=0)
	    {
	      Ngrid[j]=Ngrid[ind[i]];
	      Mgrid_l[j]=Mgrid_l[ind[i]];
	      Mgrid_u[j]=Mgrid_u[ind[i]];
	      //Rprintf("%d\t%lf\t%lf\t%lf\n",j,Ngrid[j],Mgrid_l[j],Mgrid_u[j]);
	      j++;
	    }
	}

      k=k+newG;
      //Rprintf("added\t%d\tgrids, new number of grids\t%d\n",newG,k);

    }
  
  REAL(res)[0]=k;
  for(i=0;i<k;i++)
    {
      REAL(res)[1+7*i]=0;
      REAL(res)[2+7*i]=0;
      REAL(res)[3+7*i]=0;
      REAL(res)[4+7*i]=0;
      REAL(res)[7+7*i]=Ngrid[i];
      REAL(res)[5+7*i]=0;
      REAL(res)[6+7*i]=0;
      for(j=0;j<n;j++)
	{
	  if(M[j]>Mgrid_l[i] & M[j]<=Mgrid_u[i])
	    {
		    REAL(res)[j+7*k+1]=i+1;
		    REAL(res)[1+7*i]+=X[j];
		    REAL(res)[2+7*i]+=X[j]*X[j];
		    REAL(res)[3+7*i]+=X[j]*X[j]*X[j];
		    REAL(res)[4+7*i]+=M[j];
		    if( X[j] < 1 )
			    REAL(res)[5+7*i]++;
		    if( Y[j] < 1 )
			    REAL(res)[6+7*i]++;
	    }
	}
      REAL(res)[1+7*i] /= Ngrid[ i ];
      REAL(res)[2+7*i] /= Ngrid[i];
      REAL(res)[3+7*i] /= Ngrid[i];
      REAL(res)[4+7*i] /= Ngrid[i];
    }

  UNPROTECT(1);
  return(res);

}
