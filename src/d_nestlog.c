#include "header.h"


void expfunction_nestlog(double *x, double *par, double *result)
{
  double alpha0 = par[0]; double alpha12 = par[1];

  double W = pow(x[0], -1./(alpha12*alpha0) ) + 
    pow(x[1] , -1./(alpha12*alpha0 )) ;
  double V = pow(W, alpha12) + pow(x[2], -1./alpha0 );
  
  *result=pow(V, alpha0);

}

double logd_nestlog_point(double *x, double alpha0, double alpha12 )
{
  double W = pow(x[0], -1./(alpha12*alpha0) ) + 
    pow(x[1] , -1./(alpha12*alpha0 )) ;
  double V = pow(W, alpha12) + pow(x[2], -1./alpha0 );

  double lConst = log(1.-alpha0) - log(alpha0) - log(3) ;
  
  double lA1 = (-1./alpha0 - 1.) *log(x[2]) + 
    (-1./(alpha12*alpha0) -1.) * (log(x[0]) + log(x[1]) ) ;

  double lA2 = (alpha12-2.)*log(W) + (alpha0-3.)*log(V) ;

  double  L1 = (2.-alpha0)/alpha0 * pow(W,alpha12)  ;
  double L2 = (1.-alpha12)/(alpha12*alpha0) * V ;

  double result = lConst + lA1 + lA2 + log(L1+L2) ;
  
  return result ;
}
			
void d_nestlog( double *x, int *pnx, int *pdim, 
		double *palpha0, double *palpha12,
		int *take_logs, int *return_vector,double *result) 
{
  int n ;
  double *index_x = x ;
  double *index_res = result ;
  int isOnSimplex ;
  double point_lDensity ;
  int onePointOut = 0 ;

  *index_res = 0 ;
  for( n=0 ; n < *pnx ; n++ )
    {
      isOnSimplex = is_on_simplex(index_x,  *pdim) ;
      if( ! isOnSimplex)
	{
	  (*index_res ) = (*take_logs) ? log(ZERO_BMAMEVT) : 0 ;
	  onePointOut = 1 ;
	}

      else
	{
	  point_lDensity =
	    logd_nestlog_point(index_x, *palpha0, *palpha12 );
	
	  if(*return_vector)
	    {
	      (*index_res ) = (*take_logs) ? 
		point_lDensity : exp(point_lDensity)  ;
	      /* index_res ++ ;*/
	    }
	  else
	    (*index_res) += point_lDensity ;
	    
	}

      index_x +=  *pdim ;
      if(*return_vector) index_res++ ;
    }
  
  if( ( !(*take_logs) ) && (! (*return_vector)) )
    *index_res = exp( *index_res) ;

  if (onePointOut &&  (! (*return_vector) ))
    (*index_res ) = (*take_logs) ?  log(ZERO_BMAMEVT) : 0  ;
}
