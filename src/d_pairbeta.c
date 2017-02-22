#include "header.h"


 double compute_logConst_pairbeta (double alpha, int dim)
{
  double lK =  log(2) + lgammafn(dim-2) - log(dim) - log(dim-1)+ 
    lgammafn(alpha * dim +1) - 
    (lgammafn(2*alpha +1) + lgammafn(alpha *(dim -2))) ;

  return(lK) ;
}


static double unNorm_logPairbetaFun(double alpha, double beta_ij, double xi,
				double xj, int dim)
{
  double A1 = (2*alpha - 1) * log(xi +xj);
  double A2 =  ((dim - 2) * alpha - dim +2) * log( 1 - xi - xj );
  double A3 = lgammafn(2*beta_ij) - 2*lgammafn(beta_ij);
  double A4 = (beta_ij-1) *  (log( xi) + log(xj) - 2*log(xi+xj) ) ;

  return(A1 + A2 + A3 + A4 ) ;
        
} 

 double unNorm_density_point_pairbeta(double alpha, double *beta, 
				      double *x, int dim )
{
  int count = 0 ;
  /*double npairs = choose( dim, 2);*/
  double unNorm_density = 0 ;
  int i,j ;
  double logBetaTerm;
  for( i =0; i < (dim -1) ; i++ )
    {
      for( j = i+1 ; j< dim; j++ )
	{
	  logBetaTerm= unNorm_logPairbetaFun( alpha,  beta[count], x[i],
					  x[j],  dim) ; 
	  unNorm_density += exp( logBetaTerm ) ;
	  count ++;
	}
    }
  
  return( unNorm_density ); 
}

 double *compute_dpairbeta(double alpha, double *beta, double *x,
			       int dim,  int nx, int take_logs)
{
  double *density = calloc(nx, sizeof(double));
  if (!density) return NULL;
  
  double lConst = compute_logConst_pairbeta ( alpha, dim) ;

  int n ;
  double *index = x ;
  double unNorm_point_dens ;
  int isOnSimplex ;
 for( n=0 ; n < nx ; n++ )
    {

      isOnSimplex = is_on_simplex(index,  dim) ;
      if( ! isOnSimplex)
	{
	  density[n]= take_logs ? ZERO_BMAMEVT : 0 ;
	}
      else
	{
	  unNorm_point_dens =
	    unNorm_density_point_pairbeta( alpha, beta, index, dim ) ;

	  density[n] =exp( lConst + log( unNorm_point_dens)) ;
	}
	  index +=  dim ;
    }

  return(density) ;
}


void d_pairbeta(double *x, int *pnx,  int *pdim,
		double *alpha, double *beta, 
		int *take_logs, int* return_vector, int *is_error,
		double *result)
{
  double *density ;
  if (! (density = compute_dpairbeta( *alpha, beta, x,
				    *pdim,   *pnx, *take_logs ) ))
    {
      *is_error = 1 ;
      return ;
    }

  int n ;

   if (*return_vector)
     {
       for (n = 0; n < *pnx; n++)
	 result[n] = *take_logs ? log(density[n]) : density[n];
     }
   
   else
     {
       *result = 0;
       for (n = 0; n < *pnx; n++) 
	 *result += log(density[n]);
       
       if (!(*take_logs)) 
	 *result = exp(*result);
     }

   free(density) ;
}
