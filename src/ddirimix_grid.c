#include "header.h"
 
void ddirimix_grid(double *Xgrid, double *Ygrid,
		   int *npoints, double *Mu,  int *nmu, double *wei,
		   double *nu, int *equi, double *result)
{
  double tempRes;
  double tempPoint[3];
  int i,j;
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

	  /* if(is_on_simplex(tempPoint, 3))  */
	  /*   { */
	  tempRes = 
	    ddirimix_point(Mu, nu, 
			   tempPoint, wei, 3, *nmu, 0);
	  *result =  tempRes ;
	  /*}*/
	  
	  ++result ;
	}

    }

}



void ddirimix_grid1D(double *Xgrid, 
		   int *npoints, double *Mu,  int *nmu, double *wei,
		   double *nu,  double *result)
{
  double tempRes;
  double tempPoint[2];
  int i;
  for(i = 0; i< *npoints; i++)
    {
      tempPoint[0] = Xgrid[i] ;
      tempPoint[1] = 1.-tempPoint[0];

      tempRes = 
	ddirimix_point(Mu, nu, 
		       tempPoint, wei, 2, *nmu, 0);
      *result =  tempRes ;
      ++result ;
    }

    

}
