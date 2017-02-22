#include "header.h"

void d_trinestlog_grid(double *Xgrid, double *Ygrid, 
		      int *npoints, double *palpha,
		      double *pbeta, int *equi , 
		      double *result)
{
  int i,j ;
  double tempRes ;
  double tempPoint[3] ;

  for(j=0; j< *npoints; j++)
    {
      for(i = 0; i< *npoints; i++)
	{
	  tempPoint[0] = Xgrid[i] ;
	  tempPoint[1] = Ygrid[j] ;

	  if(*equi)
	    {
	      transf_to_rect(tempPoint) ;
	    }

	  tempPoint[2] = 1.-tempPoint[0]-tempPoint[1] ;

	    if(is_on_simplex(tempPoint, 3)) 
	      {
		tempRes = 
		  d_trinestlog_point(tempPoint, *palpha,  pbeta, 0);
		*result =  tempRes ;
	      }
	  
	  ++result ;
	}

    }
}
