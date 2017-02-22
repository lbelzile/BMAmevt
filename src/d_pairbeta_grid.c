#include "header.h"

void d_pairbeta_grid (double *Xgrid, double *Ygrid, 
		      int *npoints, double *alpha,
		      double *beta, int *equi , 
		      double *result)
{
  /*  double *resultEnd;*/
  int i,j;
  double Const = exp(compute_logConst_pairbeta(*alpha, 3) ) ;
  double tempRes ;
  /* double tempGrid[2] ; */
  double tempPoint[3];

  for(j=0; j< *npoints ; j++)
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
		  unNorm_density_point_pairbeta(*alpha, beta, 
						tempPoint, 3);
		*result = Const * tempRes ;
	      }
	  
	  ++result ;
	}

    }
}
